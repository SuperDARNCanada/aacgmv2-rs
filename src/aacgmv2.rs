use crate::coeffs::interpolate_coeffs;
use crate::coords::{Cartesian, Geocentric, Geodetic, Geographic, Magnetic, Spherical};
use crate::igrf::IgrfCoeffs;
use crate::{
    AACGMv2Error, Method, Transform, KMAX, MAX_ALTITUDE, NUM_COORDS, NUM_FLAGS, POLY_ORDER,
    RADIUS_EARTH,
};
use chrono::{DateTime, Utc};
use std::f64::consts::PI;

/// AACGM_v2 Model.
///
/// Stores coefficients for coordinate conversions, valid for a given datetime.
pub struct Aacgmv2 {
    igrf: IgrfCoeffs,
    date_time: DateTime<Utc>,
    sph_harm_model: Vec<f64>,
    cint: [[[f64; NUM_FLAGS]; NUM_COORDS]; KMAX],
}
impl Aacgmv2 {
    /// Create a new model for the given datetime.
    pub fn new(dt: DateTime<Utc>) -> Result<Self, AACGMv2Error> {
        let igrf_model = IgrfCoeffs::new(dt)?;
        let sph_harm_model = interpolate_coeffs(dt)?;
        let aacgmv2 = Aacgmv2 {
            igrf: igrf_model,
            date_time: dt,
            sph_harm_model,
            cint: [[[0.0_f64; NUM_FLAGS]; NUM_COORDS]; KMAX],
        };
        Ok(aacgmv2)
    }

    /// Updates the AACGMv2 model to the given datetime.
    pub fn set_datetime(mut self, dt: DateTime<Utc>) -> Result<(), AACGMv2Error> {
        if self.date_time != dt {
            self.date_time = dt;
            self.igrf = IgrfCoeffs::new(dt)?;
            self.sph_harm_model = interpolate_coeffs(dt)?;
        }
        Ok(())
    }

    /// Retrieve the datetime that the model is valid for.
    pub fn get_datetime(&self) -> DateTime<Utc> {
        self.date_time
    }

    /// Convert to/from AACGMv2 coordinates. Returns the correct coordinates given the
    /// conversion type.
    ///
    /// All angles are given in degrees.
    ///
    /// All AACGM-v2 conversions are done in geocentric coordinates using a value of 6371.2 km for the Earth radius.
    pub fn convert(
        &mut self,
        lat: f64,
        lon: f64,
        height: f64,
        transform: &Transform,
        method: &Method,
    ) -> Result<(f64, f64, f64), AACGMv2Error> {
        let order: usize = 10; // pass in so a lower order would be allowed?

        let mut out_lat: f64;
        let out_lon: f64;
        let mut altitude: f64;
        let mut in_lat = lat;
        let mut in_lon = lon;

        if in_lat.abs() > 90. {
            return Err(AACGMv2Error::Coords(format!(
                "latitude must be in the range -90 to +90 degrees: {in_lat}"
            )));
        }

        // if forward calculation (G2A) and input coordinates are given in geodetic coordinates,
        // then must first convert to geocentric coordinates
        match transform {
            Transform::GeodeticToAACGMv2 => {
                // modify lat/lon/alt to geocentric values
                let coords = Spherical {
                    r: height,
                    t: (90.0 - in_lat).to_radians(),
                    p: in_lon.to_radians(),
                };
                let geod = Geodetic { coords };
                let geoc = geod.to_geocentric();
                in_lat = geoc.lat();
                in_lon = geoc.lon();
                altitude = geoc.alt();
            }
            _ => {
                altitude = height;
            }
        }

        if altitude < 0.0 {
            return Err(AACGMv2Error::Coords(
                "coordinate transformations are not intended for altitudes < 0 km".to_string(),
            ));
        }

        /* altitude > 2000 km not allowed for coefficients */
        if let Method::Coeffs = method {
            if altitude > MAX_ALTITUDE {
                return Err(AACGMv2Error::Coords(format!(
                    "coefficients are not valid for altitudes above {MAX_ALTITUDE} km; You must either use \
                    field-line tracing (TRACE or ALLOWTRACE) or indicate that you know this is a very \
                    bad idea (BADIDEA)"
                )));
            }
        }

        (out_lat, out_lon) =
            self.convert_geo_coord_v2(in_lat, in_lon, altitude, transform, method, order)?;

        match transform {
            Transform::GeodeticToAACGMv2 | Transform::GeocentricToAACGMv2 => {
                altitude = (altitude + RADIUS_EARTH) / RADIUS_EARTH;
            }
            Transform::AACGMv2ToGeodetic => {
                let coords = Spherical {
                    r: (RADIUS_EARTH + altitude) / RADIUS_EARTH,
                    t: (90.0 - out_lat).to_radians(),
                    p: out_lon.to_radians(),
                };
                let geoc = Geocentric { coords };
                let geod = geoc.to_geodetic();
                out_lat = geod.coords.t.to_degrees();
                altitude = geod.coords.r;
            }
            Transform::AACGMv2ToGeocentric => {}
        }

        Ok((out_lat, out_lon, altitude))
    }

    /// Select the correct function based on `convert_type`
    fn aacgm_dispatch(
        &self,
        lat: f64,
        lon: f64,
        height: f64,
        convert_type: &Transform,
    ) -> Result<(f64, f64), AACGMv2Error> {
        match convert_type {
            Transform::GeodeticToAACGMv2 | Transform::GeocentricToAACGMv2 => {
                self.aacgmv2_trace(lat, lon, height)
            }
            Transform::AACGMv2ToGeodetic | Transform::AACGMv2ToGeocentric => {
                self.aacgmv2_trace_inv(lat, lon, height)
            }
        }
    }

    /// Second-level function used to determine the lat/lon of the input coordinates.
    ///
    /// Input Arguments:
    ///     in_lat        - latitude in degrees
    ///     in_lon        - longitude in degrees
    ///     height        - height above Earth in km
    ///     args          - ConvertArgs struct
    ///     order         - integer order of spherical harmonic expansion
    ///
    /// Return Value:
    ///     (output latitude, output longitude) in degrees
    fn convert_geo_coord_v2(
        &mut self,
        lat_in: f64,
        lon_in: f64,
        height_in: f64,
        dir: &Transform,
        method: &Method,
        order: usize,
    ) -> Result<(f64, f64), AACGMv2Error> {
        let mut ylmval: [f64; KMAX] = [0.0; KMAX];

        let colat_temp: f64;

        let ztmp: f64;
        let fac: f64;

        // Using field-line tracing
        match method {
            Method::Trace => {
                let out_coords = self.aacgm_dispatch(lat_in, lon_in, height_in, dir)?;
                return Ok(out_coords);
            }
            Method::AllowTrace => {
                if height_in > MAX_ALTITUDE {
                    let out_coords = self.aacgm_dispatch(lat_in, lon_in, height_in, dir)?;
                    return Ok(out_coords);
                }
            }
            _ => {}
        }

        // Using coefficients, not field-line tracing

        // determine the altitude dependence of the coefficients
        let flag: usize = match dir {
            Transform::GeocentricToAACGMv2 | Transform::GeodeticToAACGMv2 => 0,
            Transform::AACGMv2ToGeocentric | Transform::AACGMv2ToGeodetic => 1,
        };
        let alt_var = height_in / MAX_ALTITUDE;
        let alt_var_sq = alt_var * alt_var;
        let alt_var_cu = alt_var * alt_var_sq;
        let alt_var_qu = alt_var * alt_var_cu;

        for i in 0..NUM_COORDS {
            for j in 0..KMAX {
                let offset = j + i * KMAX + flag * KMAX * NUM_COORDS * POLY_ORDER;
                let multiplier = KMAX * NUM_COORDS;
                // TODO: change to allow general polynomial approximation
                self.cint[j][i][flag] = self.sph_harm_model[offset]
                    + self.sph_harm_model[offset + multiplier] * alt_var
                    + self.sph_harm_model[offset + 2 * multiplier] * alt_var_sq
                    + self.sph_harm_model[offset + 3 * multiplier] * alt_var_cu
                    + self.sph_harm_model[offset + 4 * multiplier] * alt_var_qu;
            }
        }

        let lon_input = lon_in.to_radians();

        let colat_input = if flag == 0 {
            (90. - lat_in).to_radians()
        } else {
            // use intermediate "at-altitude" coordinates for inverse transform
            let lat_adj = cgm_to_alt(height_in, lat_in)?;
            (90. - lat_adj).to_radians()
        };

        /* Compute the values of the spherical harmonic functions.
         * NOTE: this function was adapted to use orthonormal SH functions */
        compute_sph_harmonics(colat_input, lon_input, order, &mut ylmval);

        let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);

        for (k, val) in ylmval.iter().enumerate().take((order + 1) * (order + 1)) {
            x += self.cint[k][0][flag] * val;
            y += self.cint[k][1][flag] * val;
            z += self.cint[k][2][flag] * val;
        }

        /* COMMENT: SGS
         *
         * This answers one of my questions about how the coordinates for AACGM are
         * guaranteed to be on the unit sphere. Here they compute xyz independently
         * using the SH coefficients for each coordinate. They reject anything that
         * is +/- .1 Re from the surface of the Earth. They then scale each xyz
         * coordinate by the computed radial distance. This is a TERRIBLE way to do
         * things... but necessary for the inverse transformation.
         */
        /* SGS - new method that ensures position is on unit sphere and results in a
         *       much better fit. Uses z coordinate only for sign, i.e., hemisphere.
         */
        if flag == 0 {
            fac = x * x + y * y;
            if fac > 1. {
                /* we are in the forbidden region and the solution is undefined */
                return Err(AACGMv2Error::Internal(
                    -64,
                    "Position too far from unit sphere".to_string(),
                ));
            }

            ztmp = (1. - fac).sqrt();
            z = if z < 0.0 { -ztmp } else { ztmp };

            colat_temp = z.acos();
        } else {
            /* SGS - for inverse the old normalization produces lower overall errors...*/
            let r = (x * x + y * y + z * z).sqrt();
            if !(0.9..=1.1).contains(&r) {
                return Err(AACGMv2Error::Internal(
                    -32,
                    "Position too far from unit sphere".to_string(),
                ));
            }

            z /= r;
            x /= r;
            y /= r;

            if z > 1. {
                colat_temp = 0.0;
            } else if z < -1. {
                colat_temp = PI;
            } else {
                colat_temp = z.acos();
            }
        }

        let lon_temp = if (x.abs() < 1e-8) && (y.abs() < 1e-8) {
            0.0
        } else {
            y.atan2(x)
        };

        let lat_out = 90. - colat_temp.to_degrees();
        let lon_out = lon_temp.to_degrees();

        Ok((lat_out, lon_out))
    }

    fn aacgmv2_trace(
        &self,
        lat_in: f64,
        lon_in: f64,
        alt: f64,
    ) -> Result<(f64, f64), AACGMv2Error> {
        let mut rtp = Spherical::default();
        let mut geo_coords: Geographic;
        let mut mag_coords: Magnetic;
        let mut current_coords: Geographic;
        let mut previous_coords = Geographic::default();

        // Q: these could eventually be command-line options
        let mut step_size_RE = 1.0 / RADIUS_EARTH;
        let step_size_RE_init = step_size_RE;
        let eps = 1.0e-4 / RADIUS_EARTH;

        // for the model we are doing the tracing in geocentric coordinates
        rtp.r = (RADIUS_EARTH + alt) / RADIUS_EARTH; // distance in RE; 1.0 is surface of sphere
        rtp.t = (90. - lat_in).to_radians(); // colatitude in radians
        rtp.p = lon_in.to_radians(); // longitude in radians

        // convert position to Cartesian coords
        geo_coords = rtp.to_cartesian().to_geo();

        // convert to magnetic Dipole coordinates
        mag_coords = self.igrf.geo2mag(geo_coords);

        let idir = if mag_coords.coords.z > 0. { -1 } else { 1 }; // N or S hemisphere

        step_size_RE = step_size_RE_init;

        /*
        ; trace to magnetic equator
        ;
        ; Note that there is the possibility that the magnetic equator lies
        ; at an altitude above the surface of the Earth but below the starting
        ; altitude. I am not certain of the definition of CGM, but these
        ; fieldlines map to very different locations than the solutions that
        ; lie above the starting altitude. I am considering the solution for
        ; this set of fieldlines as undefined; just like those that lie below
        ; the surface of the Earth.
        ;
        ; Added a check for when tracing goes below altitude so as not to continue
        ; tracing beyond what is necessary.
        ;
        ; Also making sure that stepsize does not go to zero
        */
        let mut below = false;
        let mut k = 0;
        while !below && (idir as f64 * mag_coords.coords.z < 0.) {
            previous_coords = geo_coords; // initialize

            (geo_coords.coords, step_size_RE) =
                self.runge_kutta_45(geo_coords.coords, idir, step_size_RE, eps, 1)?; // set to 0 for RK4

            // make sure that stepsize does not go to zero
            if step_size_RE * RADIUS_EARTH < 1e-2 {
                step_size_RE = 1e-2 / RADIUS_EARTH;
            }

            // convert to magnetic Dipole coordinates
            mag_coords = self.igrf.geo2mag(geo_coords);
            below = geo_coords.norm().powi(2)
                < (RADIUS_EARTH + alt) * (RADIUS_EARTH + alt) / (RADIUS_EARTH * RADIUS_EARTH);
            k += 1;
        }
        let niter = k;

        if !below && niter > 1 {
            // now bisect stepsize (fixed) to land on magnetic equator within 1 m
            current_coords = previous_coords;
            while step_size_RE > 1e-3 / RADIUS_EARTH {
                step_size_RE *= 0.5;
                previous_coords = current_coords;
                (current_coords.coords, step_size_RE) =
                    self.runge_kutta_45(current_coords.coords, idir, step_size_RE, eps, 0)?; // using RK4
                mag_coords = self.igrf.geo2mag(current_coords);

                // Is it possible that resetting here causes a doubling of the tol?
                if idir as f64 * mag_coords.coords.z > 0.0 {
                    current_coords = previous_coords;
                }
            }
        } else {
            current_coords = geo_coords;
        }

        // 'trace' back to reference surface along Dipole field lines
        let Lshell = current_coords.norm();
        let lat_out = -idir as f64 * (1. / Lshell).sqrt().acos().to_degrees();
        let mut lon_out: f64;
        if Lshell < (RADIUS_EARTH + alt) / RADIUS_EARTH {
            return Err(AACGMv2Error::Internal(
                -1,
                format!(
                    "Magnetic equator below the given altitude ({Lshell} < {}), alt={alt}",
                    (RADIUS_EARTH + alt) / RADIUS_EARTH
                ),
            ));
        } else {
            mag_coords = self.igrf.geo2mag(current_coords); // geographic to magnetic
            rtp = mag_coords.to_spherical();

            lon_out = rtp.p.to_degrees();
            if lon_out > 180.0 {
                lon_out -= 360.0;
            }
        }

        Ok((lat_out, lon_out))
    }

    fn aacgmv2_trace_inv(
        &self,
        mut lat_in: f64,
        lon_in: f64,
        alt: f64,
    ) -> Result<(f64, f64), AACGMv2Error> {
        let mut rtp: Spherical;
        let mut geo_coords: Geographic;
        let mut mag_coords = Magnetic::default();
        let mut current_coords: Geographic;
        let mut previous_coords = Geographic::default();

        /* Q: these could eventually be command-line options */
        let mut step_size_RE = 1.0 / RADIUS_EARTH;
        let eps = 1.0e-4 / RADIUS_EARTH;

        /* Q: Test this */
        /* poles map to infinity */
        if (lat_in.abs() - 90.).abs() < 1.0e-6 {
            let delta = if lat_in > 0.0 { -1.0e-6 } else { 1.0e-6 };
            lat_in += delta;
        }

        let (lat_out, mut lon_out): (f64, f64);
        let Lshell = 1. / (lat_in.to_radians().cos() * lat_in.to_radians().cos());
        if Lshell < (RADIUS_EARTH + alt) / RADIUS_EARTH {
            /* solution does not exist; the starting
             * position at the magnetic equator is below
             * the altitude of interest */
            return Err(AACGMv2Error::Internal(
                -1,
                "Magnetic equator below the given altitude".to_string(),
            ));
        }

        /* magnetic Cartesian coordinates of fieldline trace starting point */
        mag_coords.coords.x = Lshell * lon_in.to_radians().cos();
        mag_coords.coords.y = Lshell * lon_in.to_radians().sin();
        mag_coords.coords.z = 0.;

        /* geographic Cartesian coordinates of starting point */
        geo_coords = self.igrf.mag2geo(mag_coords);

        /* geographic spherical coordinates of starting point */
        rtp = geo_coords.coords.to_spherical();

        let mut num_iter = 0;

        // direction of trace is determined by the starting hemisphere?
        let idir = if lat_in > 0.0 { 1 } else { -1 };

        /* trace back to altitude above Earth */
        while rtp.r > (RADIUS_EARTH + alt) / RADIUS_EARTH {
            previous_coords = geo_coords;

            (geo_coords.coords, step_size_RE) =
                self.runge_kutta_45(geo_coords.coords, idir, step_size_RE, eps, 1)?; /* set to 0 for RK4: /noadapt)*/

            // make sure that stepsize does not go to zero
            if step_size_RE * RADIUS_EARTH < 5.0e-1 {
                step_size_RE = 5e-1 / RADIUS_EARTH;
            }

            rtp = geo_coords.coords.to_spherical();
            num_iter += 1;
        }

        if num_iter > 1 {
            // now bisect stepsize (fixed) to land on magnetic equator w/in 1 m
            current_coords = previous_coords;
            while step_size_RE > 1e-3 / RADIUS_EARTH {
                step_size_RE *= 0.5;
                previous_coords = current_coords;
                (current_coords.coords, step_size_RE) =
                    self.runge_kutta_45(current_coords.coords, idir, step_size_RE, eps, 0)?; // using RK4

                rtp = current_coords.coords.to_spherical();
                if rtp.r < (RADIUS_EARTH + alt) / RADIUS_EARTH {
                    current_coords = previous_coords;
                }
            }
        }

        lat_out = rtp.lat();
        lon_out = rtp.lon();
        if lon_out > 180.0 {
            lon_out -= 360.;
        }

        Ok((lat_out, lon_out))
    }

    /// Advance position along magnetic field line by one step, i.e.,
    /// numerical field-line tracing using either a fixed stepsize RK4 method
    /// or a Runge-Kutta-Fehlberg adaptive stepsize ODE solver.
    ///
    /// Input Arguments:
    ///     xyz           - Cartesian position
    ///     dir           - direction along field-line to trace
    ///     ds            - stepsize to take
    ///
    /// Keywords:
    ///     fixed         - set this keyword to do RK4 method with stepsize ds
    ///     max_ds        - maximum stepsize that is allowed, in units of RE
    ///     RRds          - set to use a maximum stepsize that is proportional
    ///                     to cube of the distance from the origin.
    fn runge_kutta_45(
        &self,
        start: Cartesian,
        idir: i32,
        mut step_size: f64,
        eps: f64,
        code: i32,
    ) -> Result<(Cartesian, f64), AACGMv2Error> {
        /* convert position to spherical coords */
        let rtp = start.to_spherical();

        /* compute IGRF field in spherical coords */
        let mag_field_sph = self.igrf.compute(&rtp)?;

        /* convert field from spherical coords to Cartesian */
        let b_cartesian = rtp.vec_to_cartesian(mag_field_sph);

        /* magnitude of field to normalize vector */
        let b_field_mag = b_cartesian.norm();

        if code == 0 {
            // no adaptive stepping
            // RK4 Method
            let k1 = step_size * idir as f64 * b_cartesian / b_field_mag;
            let temp = start + 0.5 * k1;
            let k2 = self.step_along_field_line(temp, idir, step_size)?;
            let temp = start + 0.5 * k2;
            let k3 = self.step_along_field_line(temp, idir, step_size)?;
            let temp = start + k3;
            let k4 = self.step_along_field_line(temp, idir, step_size)?;
            let end = start + (k1 + k2 + k2 + k3 + k3 + k4) / 6.0;
            return Ok((end, step_size));
        }
        /************************\
         * Adaptive RK45 method *
        \************************/
        let mut rr = eps + 1.0; /* just to get into the loop */
        let mut start_pos = start;
        while rr > eps {
            let k1 = step_size * idir as f64 * b_cartesian / b_field_mag;
            let temp = start + k1 / 4.0;
            let k2 = self.step_along_field_line(temp, idir, step_size)?;
            let temp = start + (3.0 * k1 + 9.0 * k2) / 32.0;
            let k3 = self.step_along_field_line(temp, idir, step_size)?;
            let temp = start + (1932. * k1 - 7200. * k2 + 7296. * k3) / 2197.;
            let k4 = self.step_along_field_line(temp, idir, step_size)?;
            let temp = start + 439. * k1 / 216. - 8. * k2 + 3680. * k3 / 513. - 845. * k4 / 4104.;
            let k5 = self.step_along_field_line(temp, idir, step_size)?;
            let temp = start - 8. * k1 / 27. + 2. * k2 - 3544. * k3 / 2565. + 1859. * k4 / 4104.
                - 11. * k5 / 40.;
            let k6 = self.step_along_field_line(temp, idir, step_size)?;

            let w1 = start + 25. * k1 / 216. + 1408. * k3 / 2565. + 2197. * k4 / 4104. - k5 / 5.;
            let w2 = start + 16. * k1 / 135. + 6656. * k3 / 12825. + 28561. * k4 / 56430.
                - 9. * k5 / 50.
                + 2. * k6 / 55.;

            rr = (w1 - w2).norm() / step_size;
            if rr.abs() > 1.0e-16 {
                let delt = 0.84 * (eps / rr).powf(0.25); /* this formula sucks because I have
                                                         no it where it came from.
                                                         Obviously it involves factors in
                                                         the LTEs of the two methods, but
                                                         I cannot find them written down
                                                         anywhere. */

                /*newds = ds * delt;
                //ds = newds;*/
                step_size *= delt;
                /* maximum stepsize is fixed to max_ds in units of Re */
                /*if keyword_set(max_ds) then ds = min([max_ds,ds])*/
                /* maximum stepsize is r^2 * 1km, where r is in units of Re */
                /*if keyword_set(RRds) then   ds = min([50*r*r*r/RE, ds])*/
                step_size = (50.0 * rtp.r * rtp.r * rtp.r / RADIUS_EARTH).min(step_size);
            } /* otherwise leave the stepsize alone */

            /* we use the RK4 solution */
            start_pos = w1;
            /*
            ; I would assume that using the higher order RK5 method is better, but
            ; there is the suggestion that using the RK4 solution guarantees accuracy
            ; while the RK5 does not. Apparently some texts are now suggesting using
            ; the RK5 solution...
            for (k=0;k<3;k++) xyz[k] = w2[k];
            */
        }
        let end = start_pos;

        Ok((end, step_size))
    }

    /// Advance position along magnetic field line by one step, i.e.,
    /// numerical field-line tracing.
    ///
    /// Input Arguments:
    ///     start     - Cartesian position
    ///     idir      - direction along field-line to trace
    ///     step_size - stepsize to take, in km
    ///
    /// Errors if the IGRF model cannot compute the magnetic field at `start`.
    fn step_along_field_line(
        &self,
        start: Cartesian,
        idir: i32,
        step_size: f64,
    ) -> Result<Cartesian, AACGMv2Error> {
        let sph = start.to_spherical();
        let mag_sph = self.igrf.compute(&sph)?;
        let mag_cart = sph.vec_to_cartesian(mag_sph); /* convert field to Cartesian */
        let bmag = mag_cart.norm();

        let end = (step_size * idir as f64 / bmag) * mag_cart;

        Ok(end)
    }
}

/// Transformation from AACGM to so-called 'at-altitude' coordinates.
/// The purpose of this function is to scale the latitudes in such a
/// way so that there is no gap. The problem is that for non-zero
/// altitudes (h) are range of latitudes near the equator lie on dipole
/// field lines that near reach the altitude h, and are therefore not
/// accessible. This mapping closes the gap.
///
/// Uses the formula `cos (lat_at-alt) = sqrt( (Re + h)/Re ) cos (lat_aacgm)`
///
/// Input Arguments:
///     r_height_in   - The altitude (h)
///     r_lat_in      - The AACGM latitude
///
/// Returns:
///     r_lat_adj     - The 'at-altitude' latitude
///     error         - variable is set if latitude is below the value that
///                     is mapped to the origin
fn cgm_to_alt(r_height_in: f64, r_lat_in: f64) -> Result<f64, AACGMv2Error> {
    /* convert from AACGM to at-altitude coordinates */
    let mut r1 = r_lat_in.to_radians().cos();
    let ra = (r_height_in / RADIUS_EARTH + 1.0) * (r1 * r1);
    if ra > 1.0 {
        return Err(AACGMv2Error::Internal(
            -64,
            "Unable to convert AACGM latitude to `at-altitude` latitude".to_string(),
        ));
    }

    r1 = ra.sqrt().acos();
    let r_lat_adj = if r1 * r_lat_in < 0.0 {
        -r1.to_degrees()
    } else {
        r1.to_degrees()
    };

    Ok(r_lat_adj)
}

/// Computes an array of real spherical harmonic function values
/// Y_lm(phi,theta) for a given colatitiude (phi) and longitude (theta)
/// for all the values up to l = order, which is typically 10. The
/// values are stored in a 1D array of dimension (order+1)^2. The
/// indexing scheme used is (l, -l), (l, -l+1), ... (l, l-1), (l, l), and so forth.
///
/// Input Arguments:
///     colat: The colatitude of the point for which the spherical
///            harmonic Y_lm is to be calculated
///     lon: The longitude of the point for which the spherical
///          harmonic Y_lm is to be calculated
///     order: The order of the spherical harmonic function expansion.
///            The total number of terms computed will be (order+1)^2
///
/// NOTES by SGS:
///
/// It is likely that the original version was taken from FORTRAN and used array
/// indexing that begins with 1. Indexing is somewhat more natural using the
/// zeros-based indexing of C/IDL. Indices have thus been changed from the
/// original version.
///
/// It appears that the original version used unnormalized spherical harmonic
/// functions. I suspect this might be better, but realized it too late. The
/// coefficients I derived are for orthonormal spherical harmonic functions
/// which then require the same for evaluation. I believe that the original
/// authors used orthogonal spherical harmonic functions which eliminate the
/// need for computing the normalization factors. I suspect this is just fine,
/// but have not tested it.
fn compute_sph_harmonics(colat: f64, lon: f64, order: usize, ylmval: &mut [f64]) {
    let cos_theta = colat.cos();
    let sin_theta = colat.sin();

    let cos_lon = lon.cos();
    let sin_lon = lon.sin();

    let d1 = -sin_theta;
    let mut z2_r = cos_lon;
    let mut z2_i = sin_lon;

    let mut z1_r = d1 * z2_r;
    let mut z1_i = d1 * z2_i;
    let q_fac_r = z1_r;
    let q_fac_i = z1_i;

    /*
     * Generate Zonal Harmonics (P_l^(m=0) for l = 1,order) using recursion
     * relation (6.8.7), p. 252, Numerical Recipes in C, 2nd. ed., Press. W.
     * et al. Cambridge University Press, 1992) for case where m = 0.
     *
     * l Pl = cos(theta) (2l-1) Pl-1 - (l-1) Pl-2          (6.8.7)
     *
     * where Pl = P_l^(m=0) are the associated Legendre polynomials
     *
     */

    ylmval[0] = 1.0; /* l = 0, m = 0 */
    ylmval[2] = cos_theta; /* l = 1, m = 0 */

    for l in 2..order + 1 {
        /* indices for previous two values: k = l * (l+1) + m with m=0 */
        let ia = (l - 2) * (l - 1);
        let ib = (l - 1) * l;
        let ic = l * (l + 1);

        ylmval[ic] =
            (cos_theta * (2 * l - 1) as f64 * ylmval[ib] - (l - 1) as f64 * ylmval[ia]) / l as f64;
    }

    /*
     * Generate P_l^l for l = 1 to (order+1)^2 using algorithm based upon (6.8.8)
     * in Press et al., but incorporate longitude dependence, i.e., sin/cos (phi)
     *
     * Pll = (-1)^l (2l-1)!! (sin^2(theta))^(l/2)
     *
     * where Plm = P_l^m are the associated Legendre polynomials
     *
     */

    let mut q_val_r = q_fac_r;
    let mut q_val_i = q_fac_i;
    ylmval[3] = q_val_r; /* l = 1, m = +1 */
    ylmval[1] = -q_val_i; /* l = 1, m = -1 */
    for l in 2..order + 1 {
        let d1 = (l * 2 - 1) as f64;
        z2_r = d1 * q_fac_r;
        z2_i = d1 * q_fac_i;
        z1_r = z2_r * q_val_r - z2_i * q_val_i;
        z1_i = z2_r * q_val_i + z2_i * q_val_r;
        q_val_r = z1_r;
        q_val_i = z1_i;

        /* indices for previous two values: k = l * (l+1) + m */
        let ia = l * (l + 2); /* m = +l */
        let ib = l * l; /* m = -l */

        ylmval[ia] = q_val_r;
        ylmval[ib] = -q_val_i;
    }

    /*
     * Generate P_l,l-1 to P_(order+1)^2,l-1 using algorithm based upon (6.8.9)
     * in Press et al., but incorporate longitude dependence, i.e., sin/cos (phi)
     *
     * Pl,l-1 = cos(theta) (2l-1) Pl-1,l-1
     *
     */

    for l in 2..order + 1 {
        let l2 = l * l;
        let tl = 2 * l;
        /* indices for Pl,l-1; Pl-1,l-1; Pl,-(l-1); Pl-1,-(l-1) */
        let ia = l2 - 1;
        let ib = l2 - tl + 1;
        let ic = l2 + tl - 1;
        let id = l2 + 1;

        let fac = tl - 1;
        ylmval[ic] = fac as f64 * cos_theta * ylmval[ia]; /* Pl,l-1   */
        ylmval[id] = fac as f64 * cos_theta * ylmval[ib]; /* Pl,-(l-1) */
    }

    /*
     * Generate remaining P_l+2,m to P_(order+1)^2,m for each m = 1 to order-2
     * using algorithm based upon (6.8.7) in Press et al., but incorporate
     * longitude dependence, i.e., sin/cos (phi).
     *
     * for each m value 1 to order-2 we have P_mm and P_m+1,m so we can compute
     * P_m+2,m; P_m+3,m; etc.
     *
     */

    for m in 1..order - 1 {
        for l in m + 2..order + 1 {
            let ca = ((2 * l - 1) as f64) / (l - m) as f64;
            let cb = ((l + m - 1) as f64) / (l - m) as f64;

            let l2 = l * l;
            let mut ic = l2 + l + m;
            let mut ib = l2 - l + m;
            let mut ia = l2 - l - l - l + 2 + m;
            /* positive m */
            ylmval[ic] = ca * cos_theta * ylmval[ib] - cb * ylmval[ia];

            ic -= m + m;
            ib -= m + m;
            ia -= m + m;
            /* negative m */
            ylmval[ic] = ca * cos_theta * ylmval[ib] - cb * ylmval[ia];
        }
    }

    /*
     * Normalization added here (SGS)
     *
     * Note that this is NOT the standard spherical harmonic normalization factors
     *
     * The recursive algorithms above treat positive and negative values of m in
     * the same manner. In order to use these algorithms the normalization must
     * also be modified to reflect the symmetry.
     *
     * Output values have been checked against those obtained using the internal
     * IDL legendre() function to obtain the various associated legendre
     * polynomials.
     *
     * As stated above, I think that this normalization may be unnecessary. The
     * important thing is that the various spherical harmonics are orthogonal,
     * rather than orthonormal.
     *
     */

    /* determine array of factorials */
    let mut fact = vec![0.0_f64; 2 * order + 2];

    fact[0] = 1.0;
    fact[1] = 1.0;
    for k in 2..2 * order + 2 {
        fact[k] = k as f64 * fact[k - 1];
    }

    let mut ffff = vec![0.0_f64; (order + 1) * (order + 1)];

    /* determine normalization factors */
    for l in 0..order + 1 {
        for m in 0..l + 1 {
            let k = l * (l + 1) + m; /* 1D index for l,m */
            ffff[k] = ((2 * l + 1) as f64 / (4.0 * PI) * fact[l - m] / fact[l + m]).sqrt();
            ylmval[k] *= ffff[k];
        }
        for m in -(l as i32)..0 {
            let k = (l * (l + 1)) as i32 + m; /* 1D index for l,m */
            let kk = (l * (l + 1)) as i32 - m;
            let sign = if -m % 2 == 1 { -1.0 } else { 1.0 };

            ylmval[k as usize] *= ffff[kk as usize] * sign;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use chrono::NaiveDate;

    #[test]
    fn step_along() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let model = Aacgmv2::new(dt).unwrap();

        let start = Cartesian {
            x: 0.759117255,
            y: -0.330071088,
            z: 0.837531111,
        };
        let idir = -1;
        let step_size = 1.0 / RADIUS_EARTH;
        let res = model.step_along_field_line(start, idir, step_size);
        assert!(res.is_ok());
        let end = res.unwrap();
        assert_relative_eq!(end.x, 0.000141611, max_relative = 1.0e-4);
        assert_relative_eq!(end.y, -0.000051618, max_relative = 1.0e-4);
        assert_relative_eq!(end.z, 0.000043785, max_relative = 1.0e-4);
    }

    #[test]
    fn runge_kutta() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let model = Aacgmv2::new(dt).unwrap();

        let start = Cartesian {
            x: 2.036642435,
            y: -0.662403948,
            z: 0.317608168,
        };
        let idir = 1;
        let step_size = 6.4630766e-02;
        let eps = 1.569563e-08;
        let code = 1;
        let res = model.runge_kutta_45(start, idir, step_size, eps, code);
        assert!(res.is_ok());
        let (end, step_size) = res.unwrap();

        assert_relative_eq!(end.x, 2.007697030, max_relative = 1.0e-7);
        assert_relative_eq!(end.y, -0.662095893, max_relative = 1.0e-7);
        assert_relative_eq!(end.z, 0.375375076, max_relative = 1.0e-7);
        assert_relative_eq!(step_size, 6.2884934e-02, max_relative = 1.0e-6)
    }

    #[test]
    fn geocnvrt_geod2aacgm_coeffs() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let mut model = Aacgmv2::new(dt).unwrap();

        let res = model.convert_geo_coord_v2(
            45.336703161,
            -23.5,
            1131.097495059,
            &Transform::GeodeticToAACGMv2,
            &Method::Coeffs,
            10,
        );
        assert!(res.is_ok());
        let (lat, lon) = res.unwrap();
        assert_relative_eq!(lat, 47.402896802, max_relative = 1.0e-6);
        assert_relative_eq!(lon, 56.602299929, max_relative = 1.0e-6);
    }

    #[test]
    fn geocnvrt_geod2aacgm_trace() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let mut model = Aacgmv2::new(dt).unwrap();

        let res = model.convert_geo_coord_v2(
            45.336703161,
            -23.5,
            1131.097495059,
            &Transform::GeodeticToAACGMv2,
            &Method::Trace,
            10,
        );
        assert!(res.is_ok());
        let (lat, lon) = res.unwrap();
        assert_relative_eq!(lat, 47.408677778, max_relative = 1.0e-8);
        assert_relative_eq!(lon, 56.600153826, max_relative = 1.0e-8);
    }

    #[test]
    fn geocnvrt_aacgm2geod_coeffs() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let mut model = Aacgmv2::new(dt).unwrap();

        let res = model.convert_geo_coord_v2(
            47.402896802,
            56.602299929,
            1131.097495059,
            &Transform::AACGMv2ToGeodetic,
            &Method::Coeffs,
            10,
        );
        assert!(res.is_ok());
        let (lat, lon) = res.unwrap();
        assert_relative_eq!(lat, 45.276561088, max_relative = 1.0e-6);
        assert_relative_eq!(lon, -23.477496387, max_relative = 1.0e-6);
    }

    #[test]
    fn geocnvrt_aacgm2geod_trace() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let mut model = Aacgmv2::new(dt).unwrap();

        let res = model.convert_geo_coord_v2(
            47.408677778,
            56.600153826,
            1131.097495059,
            &Transform::AACGMv2ToGeodetic,
            &Method::Trace,
            10,
        );
        assert!(res.is_ok());
        let (lat, lon) = res.unwrap();
        assert_relative_eq!(lat, 45.336703245, max_relative = 1.0e-6);
        assert_relative_eq!(lon, -23.500000111, max_relative = 1.0e-6);
    }

    #[test]
    fn test_compute_sph_harmonics() {
        let mut ylmval = [0.0_f64; KMAX];
        compute_sph_harmonics(0.824879589, 0.987896498, 10, &mut ylmval);

        assert_relative_eq!(ylmval[0], 0.282094792_f64, max_relative = 1.0e-6);
        assert_relative_eq!(ylmval[1], -0.211851367_f64, max_relative = 1.0e-6);
        assert_relative_eq!(ylmval[2], 0.331587851_f64, max_relative = 1.0e-6);
        assert_relative_eq!(ylmval[5], -0.32148387_f64, max_relative = 1.0e-6);
        assert_relative_eq!(ylmval[33], 0.424264954_f64, max_relative = 1.0e-6);
        assert_relative_eq!(ylmval[61], -0.088058059_f64, max_relative = 1.0e-6);
        assert_relative_eq!(ylmval[120], -0.022273914_f64, max_relative = 1.0e-6);
    }

    #[test]
    fn trace_inv() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let model = Aacgmv2::new(dt).unwrap();

        let res = model.aacgmv2_trace_inv(47.408678, 56.600154, 1131.0982496);
        assert!(res.is_ok());
        let (lat, lon) = res.unwrap();
        assert_relative_eq!(lat, 45.336702713, max_relative = 1.0e-6);
        assert_relative_eq!(lon, -23.499999823, max_relative = 1.0e-6);
    }
}
