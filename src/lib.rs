#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));


#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    #[test]
    fn aacgm_convert() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let mut glat: f64 = 45.5;
        let mut glon: f64 = -23.5;
        let mut height: f64 = 1135.0;
        let code = 0;

        let mut out_lat: f64 = 0.0;
        let mut out_lon: f64 = 0.0;
        let mut out_rad: f64 = 0.0;

        println!("AACGM_v2_DAT_PREFIX={:?}", std::env::var("AACGM_v2_DAT_PREFIX"));
        unsafe {
            let ret_code = crate::AACGM_v2_SetDateTime(year, month, day, hour, minute, second);
            assert_eq!(ret_code, 0);
            let ret_code = crate::AACGM_v2_Convert(glat, glon, height, &mut out_lat, &mut out_lon, &mut out_rad, code);
            assert_eq!(ret_code, 0);
        }

        assert_relative_eq!(out_lat, 47.402897, max_relative=1.0e-5);
        assert_relative_eq!(out_lon, 56.602300, max_relative=1.0e-5);
        assert_relative_eq!(out_rad, 1.177533, max_relative=1.0e-5);

        let code = 1;

        unsafe {
            let ret_code = crate::AACGM_v2_Convert(out_lat, out_lon, (out_rad - 1.0) * 6371.2, &mut glat, &mut glon, &mut height, code);
            assert_eq!(ret_code, 0);
        }

        assert_relative_eq!(glat, 45.439863, max_relative=1.0e-5);
        assert_relative_eq!(glon, -23.477496, max_relative=1.0e-5);
        assert_relative_eq!(height, 1134.977555, max_relative=1.0e-5);
    }

    #[test]
    fn aacgm_convert_trace() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let mut glat: f64 = 45.5;
        let mut glon: f64 = -23.5;
        let mut height: f64 = 1135.0;
        let code = 2;

        let mut out_lat: f64 = 0.0;
        let mut out_lon: f64 = 0.0;
        let mut out_rad: f64 = 0.0;

        println!("AACGM_v2_DAT_PREFIX={:?}", std::env::var("AACGM_v2_DAT_PREFIX"));
        unsafe {
            let ret_code = crate::AACGM_v2_SetDateTime(year, month, day, hour, minute, second);
            assert_eq!(ret_code, 0);
            let ret_code = crate::AACGM_v2_Convert(glat, glon, height, &mut out_lat, &mut out_lon, &mut out_rad, code);
            assert_eq!(ret_code, 0);
        }

        assert_relative_eq!(out_lat, 47.408678, max_relative=1.0e-5);
        assert_relative_eq!(out_lon, 56.600154, max_relative=1.0e-5);
        assert_relative_eq!(out_rad, 1.177533, max_relative=1.0e-5);

        let code = 3;

        unsafe {
            let ret_code = crate::AACGM_v2_Convert(out_lat, out_lon, (out_rad - 1.0) * 6371.2, &mut glat, &mut glon, &mut height, code);
            assert_eq!(ret_code, 0);
        }

        assert_relative_eq!(glat, 45.500000, max_relative=1.0e-5);
        assert_relative_eq!(glon, -23.500000, max_relative=1.0e-5);
        assert_relative_eq!(height, 1135.000000, max_relative=1.0e-5);
    }

    #[test]
    fn aacgm_bad_date() {
        let year = 2030;
        let month = 1;
        let day = 1;
        let hour = 0;
        let minute = 0;
        let second = 0;

        unsafe {
            let ret_code = crate::AACGM_v2_SetDateTime(year, month, day, hour, minute, second);
            assert_eq!(ret_code, -1);
        }
    }

}
