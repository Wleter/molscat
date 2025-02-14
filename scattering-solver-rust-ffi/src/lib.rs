
use std::ffi::{c_double, c_int};

use clebsch_gordan::{half_i32, half_integer::{HalfI32, HalfU32}, half_u32, wigner_3j, wigner_6j};
use quantum::states::spins::{Spin, SpinOperators};

#[no_mangle]
pub extern "C" fn hifi_mel(
    s: c_int,
    ms_c: c_int,
    ms_r: c_int,
    i: c_int,
    mi_c: c_int,
    mi_r: c_int
) -> c_double {
    let s_bra = Spin::new(
        HalfU32::from_doubled(s as u32),
        HalfI32::from_doubled(ms_c)
    );

    let s_ket = Spin::new(
        HalfU32::from_doubled(s as u32),
        HalfI32::from_doubled(ms_r)
    );

    let i_bra = Spin::new(
        HalfU32::from_doubled(i as u32),
        HalfI32::from_doubled(mi_c)
    );

    let i_ket = Spin::new(
        HalfU32::from_doubled(i as u32),
        HalfI32::from_doubled(mi_r)
    );
    
    SpinOperators::dot((s_bra,s_ket),(i_bra,i_ket))
}

#[no_mangle]
pub extern "C" fn aniso_hifi_mel(
    l: c_int,
    n_c: c_int,
    n_r: c_int,
    j_tot_c: c_int,
    j_tot_r: c_int,
    mj_tot_c: c_int,
    mj_tot_r: c_int,
    s: c_int,
    ms_c: c_int,
    ms_r: c_int,
    i: c_int,
    mi_c: c_int,
    mi_r: c_int
) -> c_double {
    let l = l as u32;

    let j_bra = n_c as u32;
    let j_ket = n_r as u32;

    let j_tot_bra = j_tot_c as u32;
    let j_tot_ket = j_tot_r as u32;

    let mr_bra = HalfI32::from_doubled(2 * mj_tot_c as i32);
    let mr_ket = HalfI32::from_doubled(2 * mj_tot_r as i32);

    let s = HalfU32::from_doubled(s as u32);
    let ms_bra = HalfI32::from_doubled(ms_c as i32);
    let ms_ket = HalfI32::from_doubled(ms_r as i32);

    let i = HalfU32::from_doubled(i as u32);
    let mi_bra = HalfI32::from_doubled(mi_c as i32);
    let mi_ket = HalfI32::from_doubled(mi_r as i32);

    let factor = (((2 * j_tot_ket + 1) * (2 * j_tot_bra + 1)
                    * (2 * j_ket + 1) * (2 * j_bra + 1)) as f64).sqrt()
                * ((2. * s.value() + 1.) * s.value() * (s.value() + 1.)
                    * (2. * i.value() + 1.) * i.value() * (i.value() + 1.)).sqrt();

    let sign = (-1.0f64).powi((j_tot_bra + j_tot_ket + l) as i32
        + ((s.double_value() + i.double_value()) as i32 
            - mr_bra.double_value() - ms_bra.double_value() - mi_bra.double_value()) / 2);

    let j_bra = HalfU32::from_doubled(2 * j_bra);
    let j_ket = HalfU32::from_doubled(2 * j_ket);
    let l = HalfU32::from_doubled(2 * l);
    let j_tot_bra = HalfU32::from_doubled(2 * j_tot_bra);
    let j_tot_ket = HalfU32::from_doubled(2 * j_tot_ket);

    let wigner = wigner_6j(j_bra, j_tot_bra, l, j_tot_ket, j_ket, half_u32!(2))
        * wigner_3j(half_u32!(1), half_u32!(1), half_u32!(2), mi_bra - mi_ket, ms_bra - ms_ket, mr_bra - mr_ket)
        * wigner_3j(j_tot_bra, half_u32!(2), j_tot_ket, -mr_bra, mr_bra - mr_ket, mr_ket)
        * wigner_3j(j_bra, half_u32!(2), j_ket, half_i32!(0), half_i32!(0), half_i32!(0))
        * wigner_3j(i, half_u32!(1), i, -mi_bra, mi_bra - mi_ket, mi_ket)
        * wigner_3j(s, half_u32!(1), s, -ms_bra, ms_bra - ms_ket, ms_ket);

    f64::sqrt(30.) / 3. * sign * factor * wigner
}

#[no_mangle]
pub extern "C" fn spin_rot_mel(
    l: c_int,
    n: c_int,
    j_tot_c: c_int,
    j_tot_r: c_int,
    mj_tot_c: c_int,
    mj_tot_r: c_int,
    s: c_int,
    ms_c: c_int,
    ms_r: c_int,
) -> c_double {
    let l = l as u32;

    let j = n as u32;

    let j_tot_bra = j_tot_c as u32;
    let j_tot_ket = j_tot_r as u32;

    let m_r_bra = HalfI32::from_doubled(2 * mj_tot_c as i32);
    let m_r_ket = HalfI32::from_doubled(2 * mj_tot_r as i32);

    let s = HalfU32::from_doubled(s as u32);
    let ms_bra = HalfI32::from_doubled(ms_c as i32);
    let ms_ket = HalfI32::from_doubled(ms_r as i32);

    let factor = (((2 * j_tot_ket + 1) * (2 * j_tot_bra + 1) 
                * (2 * j + 1) * j * (j + 1)) as f64).sqrt()
            * ((2. * s.value() + 1.) * s.value() * (s.value() + 1.)).sqrt();

    let sign = (-1.0f64).powi(1 + (j_tot_bra + j_tot_ket + l + j) as i32 - m_r_bra.double_value() / 2
                                + (s.double_value() as i32 - ms_bra.double_value()) / 2);

    let j = HalfU32::from_doubled(2 * j);
    let l = HalfU32::from_doubled(2 * l);
    let j_tot_bra = HalfU32::from_doubled(2 * j_tot_bra);
    let j_tot_ket = HalfU32::from_doubled(2 * j_tot_ket);

    let mut wigner_sum = 0.;
    for p in [half_i32!(-1), half_i32!(0), half_i32!(1)] { 
        wigner_sum += (-1.0f64).powi(p.double_value() / 2) 
            * clebsch_gordan::wigner_6j(j, j_tot_bra, l, j_tot_ket, j, half_u32!(1))
            * clebsch_gordan::wigner_3j(j_tot_bra, half_u32!(1), j_tot_ket, -m_r_bra, p, m_r_ket)
            * clebsch_gordan::wigner_3j(s, half_u32!(1), s, -ms_bra, -p, ms_ket)
    }

    factor * sign * wigner_sum
}

#[no_mangle]
pub extern "C" fn singlet_mel(
    lambda: c_int,
    l_c: c_int,
    l_r: c_int,
    n_c: c_int,
    n_r: c_int,
    j_tot: c_int,
    s: c_int,
    ms_c: c_int,
    ms_r: c_int,
    sa: c_int,
    msa_c: c_int,
    msa_r: c_int,
) -> c_double {
    let lambda = lambda as u32;
    let l_bra = l_c as u32;
    let l_ket = l_r as u32;

    let j_bra = n_c as u32;
    let j_ket = n_r as u32;

    let j_tot = j_tot as u32;

    let s = HalfU32::from_doubled(s as u32);
    let ms_bra = HalfI32::from_doubled(ms_c as i32);
    let ms_ket = HalfI32::from_doubled(ms_r as i32);

    let sa = HalfU32::from_doubled(sa as u32);
    let msa_bra = HalfI32::from_doubled(msa_c as i32);
    let msa_ket = HalfI32::from_doubled(msa_r as i32);

    let factor = (((2 * j_bra + 1) * (2 * j_ket + 1)
                * (2 * l_bra + 1) * (2 * l_ket + 1)) as f64).sqrt();

    let sign = (-1.0f64).powi((j_tot + j_bra + j_ket + s.double_value()) as i32 - sa.double_value() as i32);

    let j_bra = HalfU32::from_doubled(2 * j_bra);
    let j_ket = HalfU32::from_doubled(2 * j_ket);
    let l_bra = HalfU32::from_doubled(2 * l_bra);
    let l_ket = HalfU32::from_doubled(2 * l_ket);
    let j_tot = HalfU32::from_doubled(2 * j_tot);
    let lambda = HalfU32::from_doubled(2 * lambda);

    let wigner = clebsch_gordan::wigner_6j(j_bra, l_bra, j_tot, l_ket, j_ket, lambda)
        * clebsch_gordan::wigner_3j(j_bra, lambda, j_ket, half_i32!(0), half_i32!(0), half_i32!(0))
        * clebsch_gordan::wigner_3j(l_bra, lambda, l_ket, half_i32!(0), half_i32!(0), half_i32!(0))
        * clebsch_gordan::wigner_3j(s, sa, half_u32!(0), ms_bra, msa_bra, half_i32!(0))
        * clebsch_gordan::wigner_3j(s, sa, half_u32!(0), ms_ket, msa_ket, half_i32!(0));

    sign * factor * wigner
}

#[no_mangle]
pub extern "C" fn triplet_mel(
    lambda: c_int,
    l_c: c_int,
    l_r: c_int,
    n_c: c_int,
    n_r: c_int,
    j_tot: c_int,
    s: c_int,
    ms_c: c_int,
    ms_r: c_int,
    sa: c_int,
    msa_c: c_int,
    msa_r: c_int,
) -> c_double {
    let lambda = lambda as u32;
    let l_bra = l_c as u32;
    let l_ket = l_r as u32;

    let j_bra = n_c as u32;
    let j_ket = n_r as u32;

    let j_tot = j_tot as u32;

    let s = HalfU32::from_doubled(s as u32);
    let ms_bra = HalfI32::from_doubled(ms_c as i32);
    let ms_ket = HalfI32::from_doubled(ms_r as i32);

    let sa = HalfU32::from_doubled(sa as u32);
    let msa_bra = HalfI32::from_doubled(msa_c as i32);
    let msa_ket = HalfI32::from_doubled(msa_r as i32);

    let factor = (((2 * j_bra + 1) * (2 * j_ket + 1)
                * (2 * l_bra + 1) * (2 * l_ket + 1)) as f64).sqrt();

    let sign = (-1.0f64).powi((j_tot + j_bra + j_ket + s.double_value()) as i32 - sa.double_value() as i32);

    let j_bra = HalfU32::from_doubled(2 * j_bra);
    let j_ket = HalfU32::from_doubled(2 * j_ket);
    let l_bra = HalfU32::from_doubled(2 * l_bra);
    let l_ket = HalfU32::from_doubled(2 * l_ket);
    let j_tot = HalfU32::from_doubled(2 * j_tot);
    let lambda = HalfU32::from_doubled(2 * lambda);

    let mut triplet_wigner = 0.;
    for ms_tot in [half_i32!(-1), half_i32!(0), half_i32!(1)] {
        triplet_wigner += 3. * (-1.0f64).powi(ms_tot.double_value()) 
            * clebsch_gordan::wigner_3j(s, sa, half_u32!(1), ms_bra, msa_bra, -ms_tot)
            * clebsch_gordan::wigner_3j(s, sa, half_u32!(1), ms_ket, msa_ket, -ms_tot)
    }

    let wigner = clebsch_gordan::wigner_6j(j_bra, l_bra, j_tot, l_ket, j_ket, lambda)
        * clebsch_gordan::wigner_3j(j_bra, lambda, j_ket, half_i32!(0), half_i32!(0), half_i32!(0))
        * clebsch_gordan::wigner_3j(l_bra, lambda, l_ket, half_i32!(0), half_i32!(0), half_i32!(0));

    sign * factor * wigner * triplet_wigner
}

#[no_mangle]
pub extern "C" fn dipole_mel(
    l_c: c_int,
    l_r: c_int,
    n: c_int,
    j_tot_c: c_int,
    j_tot_r: c_int,
    mj_tot_c: c_int,
    mj_tot_r: c_int,
    s: c_int,
    ms_c: c_int,
    ms_r: c_int,
    sa: c_int,
    msa_c: c_int,
    msa_r: c_int,
) -> c_double {
    let l_bra = l_c as u32;
    let l_ket = l_r as u32;

    let n = n as u32;

    let j_tot_bra = j_tot_c as u32;
    let j_tot_ket = j_tot_r as u32;
    
    let mj_tot_bra = mj_tot_c as i32;
    let mj_tot_ket = mj_tot_r as i32;

    let s = HalfU32::from_doubled(s as u32);
    let ms_bra = HalfI32::from_doubled(ms_c as i32);
    let ms_ket = HalfI32::from_doubled(ms_r as i32);

    let sa = HalfU32::from_doubled(sa as u32);
    let msa_bra = HalfI32::from_doubled(msa_c as i32);
    let msa_ket = HalfI32::from_doubled(msa_r as i32);

    let factor = -f64::sqrt(30.) * p3(s.value()) * p3(sa.value())
        * (((2 * j_tot_bra + 1) * (2 * j_tot_ket + 1) 
            * (2 * l_bra + 1) * (2 * l_ket + 1)) as f64).sqrt();

    let sign = (-1.0f64).powi((2 * j_tot_bra + n + l_bra + l_ket) as i32 - mj_tot_bra
        + (s.double_value() as i32 - ms_bra.double_value() 
            + sa.double_value() as i32 - msa_bra.double_value()) / 2);

    let n = HalfU32::from_doubled(2 * n);
    let l_bra = HalfU32::from_doubled(2 * l_bra);
    let l_ket = HalfU32::from_doubled(2 * l_ket);
    let j_tot_bra = HalfU32::from_doubled(2 * j_tot_bra);
    let j_tot_ket = HalfU32::from_doubled(2 * j_tot_ket);
    let mj_tot_bra = HalfI32::from_doubled(2 * mj_tot_bra);
    let mj_tot_ket = HalfI32::from_doubled(2 * mj_tot_ket);

    let wigners = wigner_6j(l_bra, j_tot_bra, n, j_tot_ket, l_ket, half_u32!(2))
        * wigner_3j(j_tot_bra, half_u32!(2), j_tot_ket, -mj_tot_bra, mj_tot_bra - mj_tot_ket, mj_tot_ket)
        * wigner_3j(l_bra, half_u32!(2), l_ket, half_i32!(0), half_i32!(0), half_i32!(0))
        * wigner_3j(half_u32!(1), half_u32!(1), half_u32!(2), ms_bra - ms_ket, msa_bra - msa_ket, mj_tot_bra - mj_tot_ket)
        * wigner_3j(s, half_u32!(1), s, -ms_bra, ms_bra - ms_ket, ms_ket)
        * wigner_3j(sa, half_u32!(1), sa, -msa_bra, msa_bra - msa_ket, msa_ket);

    factor * sign * wigners
}

pub fn p3(x: f64) -> f64 {
    ((2. * x + 1.) * x * (x + 1.)).sqrt()
}

#[cfg(test)]
mod test {
    use crate::aniso_hifi_mel;

    #[test]
    fn test_aniso() {
        println!("{}", aniso_hifi_mel(1, 1, 1, 1, 1, -1, 0, 1, 1, -1, 1, 1, 1));
        println!("{}", aniso_hifi_mel(1, 1, 1, 1, 1, 0, -1, 1, -1, 1, 1, 1, 1));
    }
}