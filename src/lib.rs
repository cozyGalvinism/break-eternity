#![warn(missing_docs)]
#![crate_name = "break_eternity"]

//! A numerical library to represent numbers as large as 10^^1e308 and as 'small' as 10^-(10^^1e308).
//!
//! # Examples
//!
//! ```
//! use break_eternity::Decimal;
//!
//! let a = Decimal::from_number(1.0);
//! let b = Decimal::from_number(2.0);
//!
//! let c = a + b;
//! assert_eq!(c, Decimal::from_number(3.0));
//! ```

use custom_error::custom_error;
use std::{
    convert::{TryFrom, TryInto},
    fmt::{Display, LowerExp, UpperExp},
    num::ParseFloatError,
    ops::{Add, Div, Mul, Neg, Rem, Sub, AddAssign, SubAssign, DivAssign, MulAssign, RemAssign},
};

#[cfg(test)]
mod tests;

custom_error! {
    /// Error type for all errors in this crate.
    pub BreakEternityError
        /// An error that occurs when f_gamma and lambertw fails to converge a number (more than 100 iterations)
        IterationFailedConvering {
            /// The number that failed to converge
            z: f64
        } = "Iteration failed to converge: {z}",
        /// An error that occurs when a String cannot be parsed to a Decimal
        ParseError {
            /// The string that failed to parse
            parsed: String,
            /// The error that occurred
            error: ParseFloatError
        } = "Error while parsing \"{parsed}\": {error}",
        /// An error that occurs when lambertw is called with a number less than -1
        LambertWError = "lambertw is undefined for results < -1"
}

type Number = f64;

/// Maximum number of digits of precision to assume in a float
pub const MAX_FLOAT_PRECISION: usize = 17;
/// Exponent limit to increase a layer (9e15 is about the largest safe number)
///
/// If above this value, increase a layer
pub const EXPONENT_LIMIT: Number = 9e15;
/// Layer reduction threshold
///
/// If below this value, reduce a layer. Approximately 15.954.
pub const LAYER_REDUCTION_THRESHOLD: Number = 15.954242509439325;
/// At layer 0, smaller non-zero numbers than this become layer 1 numbers with negative mag.
///
/// After that the pattern continues as normal.
pub const FIRST_NEG_LAYER: Number = 1_f64 / 9e15;
/// Largest exponent that can appear in a number.
///
/// Not all mantissas are valid here.
pub const NUMBER_EXP_MAX: i32 = 308;
/// Smallest exponent that can appear in a number.
///
/// Not all mantissas are valid here.
pub const NUMBER_EXP_MIN: i32 = -324;
/// Amount of Es that will be printed in a string representation of a Decimal.
pub const MAX_ES_IN_A_ROW: u64 = 5;
/// Maximum number powers of 10 that will be cached.
pub const MAX_POWERS_OF_TEN: usize = (NUMBER_EXP_MAX - NUMBER_EXP_MIN + 1) as usize;
/// 2*PI
pub const TWO_PI: f64 = 6.283185307179586;
/// exp(-1)
pub const EXPN1: f64 = 0.36787944117144233;
/// W(1, 0)
pub const OMEGA: f64 = 0.5671432904097838;

lazy_static::lazy_static! {
    /// A cache of powers of 10 up to NUMBER_EXP_MAX.
    static ref POWERS_OF_TEN: Vec<f64> = {
        let mut powers_of_ten: Vec<f64> = Vec::new();
        for i in (NUMBER_EXP_MIN + 1) ..= NUMBER_EXP_MAX {
            powers_of_ten.push(format!("1e{}", i).parse().unwrap());
        }
        powers_of_ten
    };
    /// A flag if commas in a String representation should be ignored
    static ref IGNORE_COMMAS: bool = true;
    /// A flag if commas are used as decimal points instead.
    static ref COMMAS_ARE_DECIMAL_POINTS: bool = false;
}

/// Returns the sign of the number.
///
/// This implementation is different from f64::signum() because it returns 0 for 0.0 and NaN.
pub fn sign(num: f64) -> i8 {
    if num.is_nan() {
        return 0;
    }

    if num == 0.0 {
        return 0;
    }

    if num.is_infinite() {
        return if num.is_sign_positive() { 1 } else { -1 };
    }

    if num > 0.0 {
        return 1;
    }

    -1
}

/// Formats the given number to the given number of significant digits.
pub fn to_fixed(num: f64, places: i32) -> String {
    format!("{:.*}", places.try_into().unwrap(), num)
}

fn power_of_10(exp: i32) -> Number {
    POWERS_OF_TEN[(exp + 323) as usize]
}

/// Truncates the given number to the given number of significant digits and rounds, if necessary.
pub fn decimal_places(num: Number, places: i32) -> Number {
    let len = places as f64 + 1_f64;
    let num_digits = num.abs().log10().ceil();
    let rounded = (num * 10_f64.powf(len - num_digits)).round() * 10_f64.powf(num_digits - len);
    to_fixed(rounded, (len - num_digits).max(0_f64) as i32)
        .parse()
        .unwrap()
}

fn f_maglog10(num: Number) -> Number {
    sign(num) as f64 * num.abs().log10()
}

fn f_gamma(mut num: Number) -> Number {
    if !num.is_finite() {
        return num;
    }

    if num < -50.0 {
        if (num - num.trunc()).abs() < 1e-10 {
            return Number::NEG_INFINITY;
        }
        return 0.0;
    }

    let mut scal1 = 1.0;
    while num < 10.0 {
        scal1 *= num;
        num += 1.0;
    }

    num += 1.0;
    let mut l = 0.9189385332046727;
    l += (num + 0.5) * num.ln();
    l -= num;
    let num2 = num * num;
    let mut num_p = num;
    l += 1.0 / (12.0 * num_p);
    num_p *= num2;
    l += 1.0 / (360.0 * num_p);
    num_p *= num2;
    l += 1.0 / (1260.0 * num_p);
    num_p *= num2;
    l += 1.0 / (1680.0 * num_p);
    num_p *= num2;
    l += 1.0 / (1188.0 * num_p);
    num_p *= num2;
    l += 691.0 / (360360.0 * num_p);
    num_p *= num2;
    l += 7.0 / (1092.0 * num_p);
    num_p *= num2;
    l += 3617.0 / (122400.0 * num_p);

    l.exp() / scal1
}

fn f_lambertw(z: Number, mut tol: Option<Number>) -> Result<Number, BreakEternityError> {
    if tol.is_none() {
        tol = Some(1e-10);
    }

    let mut w;
    let mut wn;

    if !z.is_finite() {
        return Ok(z);
    }

    if z == 0.0 {
        return Ok(z);
    }

    if (z - 1.0).abs() < 1e-10 {
        return Ok(OMEGA);
    }

    if z < 10.0 {
        w = 0.0;
    } else {
        w = z.ln() - z.ln().ln();
    }

    for _ in 0..100 {
        wn = (z * (-w).exp() + w * w) / (w + 1.0);
        if (wn - w).abs() < tol.unwrap() * wn.abs() {
            return Ok(wn);
        } else {
            w = wn;
        }
    }

    Err(BreakEternityError::IterationFailedConvering { z })
}

fn d_lambertw(z: Decimal, mut tol: Option<Number>) -> Result<Decimal, BreakEternityError> {
    if tol.is_none() {
        tol = Some(1e-10);
    }

    let mut w;
    let mut ew;
    let mut wewz;
    let mut wn;

    if !z.mag.is_finite() {
        return Ok(z);
    }

    if z == Decimal::zero() {
        return Ok(z);
    }

    if z == Decimal::one() {
        return Ok(Decimal::from_number(OMEGA));
    }

    w = z.ln();

    // Halley's method
    for _ in 0..100 {
        ew = (-w).exp();
        wewz = w - z * ew;
        wn = w - wewz
            / (w + Decimal::from_number(1.0)
                - (w + Decimal::from_number(2.0)) * wewz
                    / (Decimal::from_number(2.0) * w + Decimal::from_number(2.0)));

        if (wn - w).abs() < Decimal::from_number(tol.unwrap()) * wn.abs() {
            return Ok(wn);
        } else {
            w = wn;
        }
    }

    Err(BreakEternityError::IterationFailedConvering { z: z.to_number() })
}

#[cfg(not(feature = "godot"))]
/// A Decimal number that can represent numbers as large as 10^^1e308 and as 'small' as 10^-(10^^1e308).
#[derive(Clone, Copy, Debug, Default)]
pub struct Decimal {
    sign: i8,
    layer: i64,
    mag: Number,
}

#[cfg(feature = "godot")]
/// A Decimal number that can represent numbers as large as 10^^1e308 and as 'small' as 10^-(10^^1e308).
#[derive(Clone, Copy, Debug, Default, gdnative::FromVariant, gdnative::ToVariant)]
pub struct Decimal {
    sign: i8,
    layer: i64,
    mag: Number,
}

impl Decimal {
    /// Returns the mantissa of the Decimal
    pub fn m(&self) -> Number {
        if self.sign == 0 {
            return 0.0;
        }

        if self.layer == 0 {
            let exp = self.mag.abs().log10().floor();
            // handle special case 5e-324
            let man = if (self.mag - 5e-324).abs() < 1e-10 {
                5.0
            } else {
                self.mag / power_of_10(exp as i32)
            };

            return self.sign as f64 * man;
        }

        if self.layer == 1 {
            let residue = self.mag - self.mag.floor();
            return self.sign as f64 * 10.0_f64.powf(residue);
        }

        self.sign as f64
    }

    /// Sets the mantissa of the Decimal
    pub fn set_m(&mut self, m: Number) {
        if self.layer <= 2 {
            self.set_from_mantissa_exponent(m, self.layer as f64);
        } else {
            // this mantissa isn't meaningful
            self.sign = sign(m);
            if self.sign == 0 {
                self.layer = 0;
                self.set_e(0.0);
            }
        }
    }

    /// Returns the exponent of the Decimal
    pub fn e(&self) -> Number {
        if self.sign == 0 {
            return 0.0;
        }

        if self.layer == 0 {
            return self.mag.log10().floor();
        }

        if self.layer == 1 {
            return self.mag.floor();
        }

        if self.layer == 2 {
            return (sign(self.mag) as f64 * self.mag.abs().powf(10.0)).floor();
        }

        self.mag * Number::INFINITY
    }

    /// Sets the exponent of the Decimal
    pub fn set_e(&mut self, e: Number) {
        self.set_from_mantissa_exponent(self.m(), e);
    }

    /// Returns the sign of the Decimal
    pub fn s(&self) -> i8 {
        self.sign
    }

    /// Sets the sign of the Decimal
    pub fn set_s(&mut self, s: i8) {
        if s == 0 {
            self.sign = 0;
            self.layer = 0;
            self.mag = 0.0;
        } else {
            self.sign = s;
        }
    }

    /// Returns the mantissa of the Decimal
    pub fn mantissa(&self) -> Number {
        self.m()
    }

    /// Sets the mantissa of the Decimal
    pub fn set_mantissa(&mut self, m: Number) {
        self.set_m(m);
    }

    /// Returns the exponent of the Decimal
    pub fn exponent(&self) -> Number {
        self.e()
    }

    /// Sets the exponent of the Decimal
    pub fn set_exponent(&mut self, e: Number) {
        self.set_e(e);
    }

    /// Creates a Decimal from a sign, a layer and a magnitude
    ///
    /// This function normalizes the inputs
    pub fn from_components(sign: i8, layer: i64, mag: Number) -> Decimal {
        Decimal::default().set_from_components(sign, layer, mag)
    }

    /// Creates a Decimal from a sign, a layer and a magnitude
    ///
    /// This function does not normalize the inputs
    pub fn from_components_no_normalize(sign: i8, layer: i64, mag: Number) -> Decimal {
        Decimal::default().set_from_components_no_normalize(sign, layer, mag)
    }

    /// Creates a Decimal from a mantissa and an exponent
    ///
    /// This function normalizes the inputs
    pub fn from_mantissa_exponent(m: Number, e: Number) -> Decimal {
        Decimal::default().set_from_mantissa_exponent(m, e)
    }

    /// Creates a Decimal from a mantissa and an exponent
    ///
    /// This function does not normalize the inputs
    pub fn from_mantissa_exponent_no_normalize(m: Number, e: Number) -> Decimal {
        Decimal::default().set_from_mantissa_exponent_no_normalize(m, e)
    }

    /// Creates a Decimal from a number (f64)
    pub fn from_number(n: Number) -> Decimal {
        if n.is_nan() {
            return Decimal::nan();
        }

        if n.is_infinite() && n.is_sign_positive() {
            return Decimal::inf();
        }

        if n.is_infinite() && n.is_sign_negative() {
            return Decimal::neg_inf();
        }

        Decimal::default().set_from_number(n)
    }

    /// Normalizes the Decimal as follows:
    ///
    /// * Whenever we are partially 0 (sign is 0 or mag and layer is 0), make it fully 0.
    /// * Whenever we are at or hit layer 0, extract sign from negative mag.
    /// * If layer === 0 and mag < FIRST_NEG_LAYER (1/9e15), shift to 'first negative layer' (add layer, log10 mag).
    /// * While abs(mag) > EXP_LIMIT (9e15), layer += 1, mag = maglog10(mag).
    /// * While abs(mag) < LAYER_DOWN (15.954) and layer > 0, layer -= 1, mag = pow(10, mag).
    /// * When we're done, all of the following should be true OR one of the numbers is not IsFinite OR layer is not IsInteger (error state):
    ///     * Any 0 is totally zero (0, 0, 0).
    ///     * Anything layer 0 has mag 0 OR mag > 1/9e15 and < 9e15.
    ///     * Anything layer 1 or higher has abs(mag) >= 15.954 and < 9e15.
    /// * We will assume in calculations that all Decimals are either erroneous or satisfy these criteria. (Otherwise: Garbage in, garbage out.)
    pub fn normalize(&mut self) -> Decimal {
        // PSEUDOCODE:
        // Whenever we are partially 0 (sign is 0 or mag and layer is 0), make it fully 0.
        // Whenever we are at or hit layer 0, extract sign from negative mag.
        // If layer === 0 and mag < FIRST_NEG_LAYER (1/9e15), shift to 'first negative layer' (add layer, log10 mag).
        // While abs(mag) > EXP_LIMIT (9e15), layer += 1, mag = maglog10(mag).
        // While abs(mag) < LAYER_DOWN (15.954) and layer > 0, layer -= 1, mag = pow(10, mag).
        // When we're done, all of the following should be true OR one of the numbers is not IsFinite OR layer is not IsInteger (error state):
        // Any 0 is totally zero (0, 0, 0).
        // Anything layer 0 has mag 0 OR mag > 1/9e15 and < 9e15.
        // Anything layer 1 or higher has abs(mag) >= 15.954 and < 9e15.
        // We will assume in calculations that all Decimals are either erroneous or satisfy these criteria. (Otherwise: Garbage in, garbage out.)

        if self.sign == 0 || (self.mag == 0.0 && self.layer == 0) {
            self.sign = 0;
            self.layer = 0;
            self.mag = 0.0;
            return *self;
        }

        if self.layer == 0 && self.mag < 0.0 {
            self.mag = -self.mag;
            self.sign = -self.sign;
        }

        if self.layer == 0 && self.mag < FIRST_NEG_LAYER {
            self.layer += 1;
            self.mag = self.mag.log10();
            return *self;
        }

        let mut abs_mag = self.mag.abs();
        let mut sign_mag = sign(self.mag) as f64;

        if abs_mag >= EXPONENT_LIMIT {
            self.layer += 1;
            self.mag = sign_mag * abs_mag.log10();
            return *self;
        }

        while abs_mag < LAYER_REDUCTION_THRESHOLD && self.layer > 0 {
            self.layer -= 1;
            if self.layer == 0 {
                self.mag = 10.0_f64.powf(self.mag);
            } else {
                self.mag = sign_mag * 10.0_f64.powf(abs_mag);
                abs_mag = self.mag.abs();
                sign_mag = sign(self.mag) as f64;
            }
        }

        if self.layer == 0 {
            if self.mag < 0.0 {
                self.mag = -self.mag;
                self.sign = -self.sign;
            } else if self.mag == 0.0 {
                self.sign = 0;
            }
        }

        *self
    }

    /// Sets the components of the Decimal from a sign, a layer and a magnitude
    ///
    /// This function normalizes the inputs
    pub fn set_from_components(&mut self, sign: i8, layer: i64, mag: Number) -> Decimal {
        self.sign = sign;
        self.layer = layer;
        self.mag = mag;

        self.normalize();
        *self
    }

    /// Sets the components of the Decimal from a sign, a layer and a magnitude
    ///
    /// This function does not normalize the inputs
    pub fn set_from_components_no_normalize(
        &mut self,
        sign: i8,
        layer: i64,
        mag: Number,
    ) -> Decimal {
        self.sign = sign;
        self.layer = layer;
        self.mag = mag;

        *self
    }

    /// Sets the components of the Decimal from a mantissa and an exponent
    ///
    /// This function normalizes the inputs
    pub fn set_from_mantissa_exponent(&mut self, m: Number, e: Number) -> Decimal {
        self.layer = 1;
        self.sign = sign(m);
        let mant = m.abs();
        self.mag = e + mant.log10();

        self.normalize();
        *self
    }

    /// Sets the components of the Decimal from a mantissa and an exponent
    ///
    /// This function does not normalize the inputs
    pub fn set_from_mantissa_exponent_no_normalize(&mut self, m: Number, e: Number) -> Decimal {
        self.set_from_mantissa_exponent(m, e);
        *self
    }

    /// Sets the components of the Decimal from a number (f64)
    pub fn set_from_number(&mut self, n: Number) -> Decimal {
        self.sign = sign(n);
        self.layer = 0;
        self.mag = n.abs();

        self.normalize();
        *self
    }

    /// Returns the Decimal as a number (f64)
    pub fn to_number(&self) -> Number {
        if self.layer == 0 {
            return self.sign as f64 * self.mag;
        }

        if self.layer == 1 {
            return self.sign as f64 * 10.0_f64.powf(self.mag);
        }

        if self.mag > 0.0 {
            if self.sign > 0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            }
        } else {
            0.0
        }
    }

    /// Returns the mantissa with the specified amount of decimal places
    pub fn mantissa_with_decimal_places(&self, places: i32) -> Number {
        if self.m().is_nan() {
            return f64::NAN;
        }

        if self.m() == 0.0 {
            return 0.0;
        }

        decimal_places(self.m(), places)
    }

    /// Returns the magnitude with the specified amount of decimal places
    pub fn magnitude_with_decimal_places(&self, places: i32) -> Number {
        if self.mag.is_nan() {
            return f64::NAN;
        }

        if self.mag == 0.0 {
            return 0.0;
        }

        decimal_places(self.mag, places)
    }

    /// Returns the Decimal as a String with the specified amount of decimal places
    pub fn to_fixed(&self, places: usize) -> String {
        if self.layer == 0 {
            return format!("{:.*}", places, self.mag);
        }

        self.to_string_with_decimal_places(places, None)
    }

    /// Returns the Decimal as a String with the specified amount of precision.
    ///
    /// If the amount of decimal places is larger than that the exponent, this will return a fixed representation of the number.
    ///
    /// Otherwise, this will return a scientific representation of the number.
    pub fn to_precision(&self, places: usize) -> String {
        if self.e() <= -7.0 {
            return format!("{:.*e}", places - 1, self);
        }

        if places as f64 > self.e() {
            return self.to_fixed(places - self.exponent() as usize - 1);
        }

        format!("{:.*e}", places - 1, self)
    }

    /// Returns the Decimal as a String with the specified amount of precision.
    ///
    /// This follows more rules than to_precision.
    ///
    /// If the layer of the Decimal is 0 and the magnitude is less than 1e21 and larger than 1e-7 (or when the magnitude is 0),
    /// this will return a fixed representation of the number.
    /// If the layer is 0, a scientific representation of the number will be returned.
    ///
    /// If the layer is 1, a scientific representation of the number will be returned.
    ///
    /// Otherwise, a scientific representation with multiple es will be returned.
    pub fn to_string_with_decimal_places(&self, places: usize, e_lower: Option<bool>) -> String {
        let e = if e_lower.unwrap_or(true) { "e" } else { "E" };

        if self.layer == 0 {
            if (self.mag < 1e21 && self.mag > 1e-7) || self.mag == 0.0 {
                return format!("{:.*}", places, self.sign as f64 * self.mag);
            }
            return format!(
                "{:.*}{}{:.*}",
                places,
                decimal_places(self.m(), places as i32),
                e,
                places,
                decimal_places(self.e(), places as i32)
            );
        }

        if self.layer == 1 {
            return format!(
                "{:.*}{}{:.*}",
                places,
                decimal_places(self.m(), places as i32),
                e,
                places,
                decimal_places(self.e(), places as i32)
            );
        }

        if self.layer <= MAX_ES_IN_A_ROW as i64 {
            return format!(
                "{}{}{:.*}",
                if self.sign > 0 { "" } else { "-" },
                e.repeat(self.layer as usize),
                places,
                decimal_places(self.mag, places as i32)
            );
        } else {
            return format!(
                "{}({}^{}){:.*}",
                if self.sign > 0 { "" } else { "-" },
                e,
                self.layer,
                places,
                decimal_places(self.mag, places as i32)
            );
        }
    }

    /// Returns the absolute value of the Decimal
    pub fn abs(&self) -> Decimal {
        Decimal::from_components_no_normalize(
            if self.sign == 0 { 0 } else { 1 },
            self.layer,
            self.mag,
        )
    }

    /// Returns a zero Decimal
    pub fn zero() -> Decimal {
        Decimal::from_components_no_normalize(0, 0, 0.0)
    }

    /// Returns a one Decimal
    pub fn one() -> Decimal {
        Decimal::from_components_no_normalize(1, 0, 1.0)
    }

    /// Returns a negative one Decimal
    pub fn neg_one() -> Decimal {
        Decimal::from_components_no_normalize(-1, 0, 1.0)
    }

    /// Returns a two Decimal
    pub fn two() -> Decimal {
        Decimal::from_components_no_normalize(1, 0, 2.0)
    }

    /// Returns a ten Decimal
    pub fn ten() -> Decimal {
        Decimal::from_components_no_normalize(1, 0, 10.0)
    }

    /// Returns a NaN Decimal
    pub fn nan() -> Decimal {
        Decimal::from_components_no_normalize(0, 0, f64::NAN)
    }

    /// Returns a positive infinity Decimal
    pub fn inf() -> Decimal {
        Decimal::from_components_no_normalize(1, 0, f64::INFINITY)
    }

    /// Returns a negative infinity Decimal
    pub fn neg_inf() -> Decimal {
        Decimal::from_components_no_normalize(-1, 0, f64::NEG_INFINITY)
    }

    /// Returns the largest safe Decimal that can be represented from an f64
    pub fn maximum() -> Decimal {
        Decimal::from_components_no_normalize(1, 0, f64::MAX)
    }

    /// Returns the smallest safe Decimal that can be represented from an f64
    pub fn minimum() -> Decimal {
        Decimal::from_components_no_normalize(1, 0, f64::MIN)
    }

    /// Rounds the Decimal to the nearest integer
    pub fn round(&self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::zero();
        }

        if self.layer == 0 {
            return Decimal::from_components(self.sign, 0, self.mag.round());
        }

        *self
    }

    /// Returns the largest Decimal less than or equal to the Decimal.
    pub fn floor(&self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::zero();
        }

        if self.layer == 0 {
            return Decimal::from_components(self.sign, 0, self.mag.floor());
        }

        *self
    }

    /// Returns the smallest Decimal greater than or equal to the Decimal.
    pub fn ceil(&self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::zero();
        }

        if self.layer == 0 {
            return Decimal::from_components(self.sign, 0, self.mag.ceil());
        }

        *self
    }

    /// Returns the integer part of the Decimal
    pub fn trunc(&self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::zero();
        }

        if self.layer == 0 {
            return Decimal::from_components(self.sign, 0, self.mag.trunc());
        }

        *self
    }

    /// Compares the absolute value of the Decimal to the absolute value of the other Decimal
    pub fn cmpabs(&self, rhs: &Decimal) -> i8 {
        let layer_a = if self.mag > 0.0 {
            self.layer as i64
        } else {
            -(self.layer as i64)
        };
        let layer_b = if rhs.mag > 0.0 {
            rhs.layer as i64
        } else {
            -(rhs.layer as i64)
        };

        if layer_a > layer_b {
            return 1;
        }

        if layer_a < layer_b {
            return -1;
        }

        if self.mag > rhs.mag {
            return 1;
        }

        if self.mag < rhs.mag {
            return -1;
        }

        0
    }

    /// Compares the absolute value of the Decimal to the absolute value of the other Decimal
    /// and returns the bigger one
    pub fn maxabs(&self, rhs: Decimal) -> Decimal {
        if self.cmpabs(&rhs) > 0 {
            *self
        } else {
            rhs
        }
    }

    /// Compares the absolute value of the Decimal to the absolute value of the other Decimal
    /// and returns the smaller one
    pub fn minabs(&self, rhs: Decimal) -> Decimal {
        if self.cmpabs(&rhs) > 0 {
            rhs
        } else {
            *self
        }
    }

    /// Returns the reciprocal of the Decimal
    pub fn recip(&self) -> Decimal {
        if self.mag == 0.0 {
            return Decimal::nan();
        }

        if self.layer == 0 {
            return Decimal::from_components(self.sign, 0, 1.0 / self.mag);
        }

        Decimal::from_components(self.sign, self.layer, -self.mag)
    }

    /// Returns the bigger of the two Decimals
    pub fn max(&self, other: Decimal) -> Decimal {
        if self > &other {
            *self
        } else {
            other
        }
    }

    /// Returns the smaller of the two Decimals
    pub fn min(&self, other: Decimal) -> Decimal {
        if self < &other {
            *self
        } else {
            other
        }
    }

    /// Clamps the Decimal to the given range
    pub fn clamp(&self, min: Decimal, max: Decimal) -> Decimal {
        self.max(min).min(max)
    }

    /// Clamps the Decimal to a minimum value
    pub fn clamp_min(&self, min: Decimal) -> Decimal {
        self.max(min)
    }

    /// Clamps the Decimal to a maximum value
    pub fn clamp_max(&self, max: Decimal) -> Decimal {
        self.min(max)
    }

    /// Tolerance is a relative tolerance, multiplied by the greater of the magnitudes of the two arguments.
    /// For example, if you put in 1e-9, then any number closer to the
    /// larger number than (larger number)*1e-9 will count as equal.
    ///
    /// Default tolerance is 1e-7
    pub fn eq_tolerance(&self, other: &Decimal, tolerance: f64) -> bool {
        if self.sign != other.sign {
            return false;
        }

        if (self.layer - other.layer).abs() > 1 {
            return false;
        }

        let mut mag_a = self.mag;
        let mut mag_b = other.mag;
        if self.layer > other.layer {
            mag_b = f_maglog10(mag_b);
        }
        if other.layer > self.layer {
            mag_a = f_maglog10(mag_a);
        }

        (mag_a - mag_b).abs() <= tolerance * mag_a.abs().max(mag_b.abs())
    }

    /// Returns the absolute log10 of the Decimal
    pub fn abs_log10(&self) -> Decimal {
        if self.sign == 0 {
            return Decimal::nan();
        }

        if self.layer > 0 {
            return Decimal::from_components(sign(self.mag), self.layer - 1, self.mag.abs());
        }

        Decimal::from_components(1, 0, self.mag.log10())
    }

    /// Returns log10 of the Decimal
    pub fn log10(&self) -> Decimal {
        if self.sign <= 0 {
            return Decimal::nan();
        }

        if self.layer > 0 {
            return Decimal::from_components(sign(self.mag), self.layer - 1, self.mag.abs());
        }

        Decimal::from_components(self.sign, 0, self.mag.log10())
    }

    /// Returns the log of the Decimal with the given base
    pub fn log(&self, base: Decimal) -> Decimal {
        if self.sign <= 0 {
            return Decimal::nan();
        }

        if base.sign <= 0 {
            return Decimal::nan();
        }

        if base.sign == 1 && base.layer == 0 && (base.mag - 1.0).abs() < 1e-10 {
            return Decimal::nan();
        }

        if self.layer == 0 && base.layer == 0 {
            return Decimal::from_components(self.sign, 0, self.mag.ln() / base.mag.ln());
        }

        self.log10() / base.log10()
    }

    /// Returns the log2 of the Decimal
    pub fn log2(&self) -> Decimal {
        if self.sign <= 0 {
            return Decimal::nan();
        }

        if self.layer == 0 {
            return Decimal::from_components(self.sign, 0, self.mag.log2());
        }

        if self.layer == 1 {
            return Decimal::from_components(
                sign(self.mag),
                0,
                self.mag.abs() * std::f64::consts::LOG2_10,
            );
        }

        if self.layer == 2 {
            return Decimal::from_components(
                sign(self.mag),
                1,
                self.mag.abs() + 0.5213902276543247,
            );
        }

        Decimal::from_components(sign(self.mag), self.layer - 1, self.mag.abs())
    }

    /// Returns the natural log of the Decimal
    pub fn ln(&self) -> Decimal {
        if self.sign <= 0 {
            return Decimal::nan();
        }

        if self.layer == 0 {
            return Decimal::from_components(self.sign, 0, self.mag.ln());
        }

        if self.layer == 1 {
            return Decimal::from_components(
                sign(self.mag),
                0,
                self.mag.abs() * std::f64::consts::LN_10,
            );
        }

        if self.layer == 2 {
            return Decimal::from_components(
                sign(self.mag),
                1,
                self.mag.abs() + 0.36221568869946325,
            );
        }

        Decimal::from_components(sign(self.mag), self.layer - 1, self.mag.abs())
    }

    /// Returns the Decimal to the power of the given exponent
    pub fn pow(self, exp: Decimal) -> Decimal {
        let a = self;
        let b = exp;

        if a.sign == 0 {
            return if b == Decimal::from_number(0.0) {
                Decimal::one()
            } else {
                a
            };
        }

        if a.sign == 1 && a.layer == 0 && (a.mag - 1.0).abs() < 1e-10 {
            return a;
        }

        if b.sign == 0 {
            return Decimal::one();
        }

        if b.sign == 1 && b.layer == 0 && (b.mag - 1.0).abs() < 1e-10 {
            return a;
        }

        let result = (a.abs_log10() * b).pow10();

        if self.sign == -1 && ((b.to_number() % 2.0).abs() - 1.0).abs() < 1e-10 {
            return -result;
        }

        result
    }

    /// Returns the Decimal raised to the next power of 10
    pub fn pow10(self) -> Decimal {
        if !self.mag.is_finite() {
            return Decimal::nan();
        }

        let mut a = self;

        if a.layer == 0 {
            let new_mag = 10.0_f64.powf(a.sign as f64 * a.mag);
            if new_mag.is_finite() && new_mag.abs() > 0.1 {
                return Decimal::from_components(1, 0, new_mag);
            } else {
                if a.sign == 0 {
                    return Decimal::one();
                }
                a = Decimal::from_components_no_normalize(a.sign, a.layer + 1, a.mag.log10());
            }
        }

        if a.sign > 0 && a.mag > 0.0 {
            return Decimal::from_components(a.sign, a.layer + 1, a.mag);
        }

        if a.sign < 0 && a.mag > 0.0 {
            return Decimal::from_components(-a.sign, a.layer + 1, -a.mag);
        }

        Decimal::one()
    }

    /// Returns the n-th root of the Decimal
    pub fn root(self, n: Decimal) -> Decimal {
        self.pow(n.recip())
    }

    /// Returns the exponential function of the Decimal
    pub fn exp(self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::one();
        }

        if self.layer == 0 && self.mag <= 709.7 {
            return Decimal::from_number((self.sign as f64 * self.mag).exp());
        }

        if self.layer == 0 {
            return Decimal::from_components(
                1,
                1,
                self.sign as f64 * std::f64::consts::E.log10() * self.mag,
            );
        }

        if self.layer == 1 {
            return Decimal::from_components(
                1,
                2,
                self.sign as f64 * std::f64::consts::LOG10_E.log10() + self.mag,
            );
        }

        Decimal::from_components(1, self.layer + 1, self.sign as f64 * self.mag)
    }

    /// Returns the gamma function of the Decimal
    pub fn gamma(&self) -> Decimal {
        if self.mag < 0.0 {
            return self.recip();
        }

        if self.layer == 0 {
            if self < &Decimal::from_components_no_normalize(1, 0, 24.0) {
                return Decimal::from_number(f_gamma(self.sign as f64 * self.mag));
            }

            let t = self.mag - 1.0;
            let mut l = 0.9189385332046727;
            l += (t + 0.5) * t.ln();
            l -= t;
            let n2 = t * t;
            let mut np = t;
            let mut lm = 12.0 * np;
            let adj = 1.0 / lm;
            let l2 = l + adj;
            if (l2 - l).abs() < 1e-10 {
                return Decimal::from_number(l).exp();
            }

            l = l2;
            np *= n2;
            lm = 1260.0 * np;
            let mut lt = 1.0 / lm;
            l += lt;
            np *= n2;
            lm = 1680.0 * np;
            lt = 1.0 / lm;
            l -= lt;
            return Decimal::from_number(l).exp();
        }

        if self.layer == 1 {
            return (*self * (self.ln() - Decimal::from_number(1.0))).exp();
        }

        self.exp()
    }

    /// Returns the factorial of the Decimal
    pub fn factorial(&self) -> Decimal {
        if self.mag < 0.0 {
            return (*self + Decimal::from_number(1.0)).gamma();
        }

        if self.layer == 0 {
            return (*self + Decimal::from_number(1.0)).gamma();
        }

        if self.layer == 1 {
            return (*self * self.ln() - Decimal::from_number(1.0)).exp();
        }

        self.exp()
    }

    /// Returns the natural logarithm of the gamma function of the Decimal
    pub fn ln_gamma(&self) -> Decimal {
        self.gamma().ln()
    }

    /// Returns the Decimal squared
    pub fn sqr(&self) -> Decimal {
        self.pow(Decimal::from_number(2.0))
    }

    /// Returns the square root of the Decimal
    pub fn sqrt(&self) -> Decimal {
        if self.layer == 0 {
            return Decimal::from_number((self.sign as f64 * self.mag).sqrt());
        }

        if self.layer == 1 {
            return Decimal::from_components(1, 2, self.mag.log10() - std::f64::consts::LOG10_2);
        }

        let mut result = Decimal::from_components_no_normalize(self.sign, self.layer - 1, self.mag)
            / Decimal::from_components_no_normalize(1, 0, 2.0);
        result.layer += 1;
        result.normalize();

        result
    }

    /// Returns the Decimal cubed
    pub fn cube(&self) -> Decimal {
        self.pow(Decimal::from_number(3.0))
    }

    /// Returns the cube root of the Decimal
    pub fn cbrt(&self) -> Decimal {
        self.pow(Decimal::from_number(1.0) / Decimal::from_number(3.0))
    }

    /// Tetrates the Decimal to the given height.
    ///
    /// Source: <https://andydude.github.io/tetration/archives/tetration2/ident.html>
    pub fn tetrate(&self, height: Option<Number>, payload: Option<Decimal>) -> Decimal {
        let mut height = height.unwrap_or(2.0_f64);
        let mut payload =
            payload.unwrap_or_else(|| Decimal::from_components_no_normalize(1, 0, 1.0));

        if height.is_infinite() && height.is_sign_positive() {
            let neg_ln = self.ln().neg();
            return neg_ln.lambertw().expect("Expected number higher than -1") / neg_ln;
        }

        if height < 0.0 {
            return payload.iteratedlog(*self, -height);
        }

        let old_height = height;
        height = height.trunc();
        let fract_height = old_height - height;

        if fract_height != 0.0 {
            if payload == Decimal::one() {
                height += 1.0;
                payload = Decimal::from_number(fract_height);
            } else if *self == Decimal::from_number(10.0) {
                payload = payload.layer_add_10(Decimal::from_number(fract_height));
            } else {
                payload = payload.layer_add(fract_height, *self);
            }
        }

        for i in 0..height as i64 {
            payload = self.pow(payload);
            // bail if we're NaN
            if !payload.mag.is_finite() {
                return payload;
            }

            if payload.layer - self.layer > 3 {
                return Decimal::from_components_no_normalize(
                    payload.sign,
                    payload.layer + (height as i64 - i - 1),
                    payload.mag,
                );
            }

            if i > 100 {
                return payload;
            }
        }

        payload
    }

    /// Returns the Decimal, iteratively exponentiated
    ///
    /// Equates to tetrating to the same height.
    pub fn iteratedexp(&self, height: Option<Number>, payload: Option<Decimal>) -> Decimal {
        self.tetrate(height, payload)
    }

    /// Iterated log: The result of applying log(base) 'times' times in a row.
    /// Approximately equal to subtratcting (times) from the number's slo representation.
    /// Equates to tetrating to a negative height.
    pub fn iteratedlog(&self, base: Decimal, mut times: f64) -> Decimal {
        if times < 0.0 {
            return base.tetrate(Some(-times), Some(*self));
        }

        let mut result = *self;
        let full_times = times;
        times = times.trunc();
        let fraction = full_times - times;

        if result.layer - base.layer > 3 {
            let layer_loss = times.min((result.layer - base.layer - 3) as f64);
            times -= layer_loss;
            result.layer -= layer_loss as i64;
        }

        for i in 0..times as i64 {
            result = result.log(base);
            if !result.mag.is_finite() {
                return result;
            }
            if i > 100 {
                return result;
            }
        }

        if fraction > 0.0 && fraction < 1.0 {
            if base == Decimal::from_number(10.0) {
                result = result.layer_add_10(Decimal::from_number(-fraction));
            } else {
                result = result.layer_add(-fraction, base);
            }
        }

        result
    }

    /// Returns the super-logarithm of the Decimal
    pub fn slog(&self, base_opt: Option<Decimal>) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::neg_one();
        }

        let mut result: f64 = 0.0;
        let base = base_opt.unwrap_or_else(|| Decimal::from_number(10.0));
        let mut copy = *self;

        if copy.layer - base.layer > 3 {
            let layer_loss = copy.layer - base.layer - 3;
            result += layer_loss as f64;
            copy.layer -= layer_loss;
        }

        for _ in 0..100 {
            if copy < Decimal::zero() {
                copy = base.pow(copy);
                result -= 1.0;
            }

            if copy <= Decimal::one() {
                return Decimal::from_number(result + copy.to_number() - 1.0);
            }

            result += 1.0;
            copy = copy.log(base);
        }

        Decimal::from_number(result)
    }

    /// Adds or removes layers from a Decimal using linear approximation
    pub fn layer_add_10(&self, diff: Decimal) -> Decimal {
        let mut diff = diff.to_number();
        let mut result = *self;

        if diff >= 1.0 {
            let layer_add = diff.trunc();
            diff -= layer_add;
            result.layer += layer_add as i64;
        }

        if diff <= -1.0 {
            let layer_add = diff.trunc();
            diff -= layer_add;
            result.layer += layer_add as i64;
            if result.layer < 0 {
                for _ in 0..100 {
                    result.layer += 1;
                    result.mag = result.mag.log10();
                    if !result.mag.is_finite() {
                        return result;
                    }

                    if result.layer >= 0 {
                        break;
                    }
                }
            }
        }

        if diff > 0.0 {
            let mut subtract_layers_later: i64 = 0;
            while result.mag.is_finite() && result.mag < 10.0 {
                result.mag = 10.0_f64.powf(result.mag);
                subtract_layers_later += 1;
            }

            if result.mag > 1e10_f64 {
                result.mag = result.mag.log10();
                result.layer += 1;
            }

            let diff_to_next_slog = (1e10_f64.ln() / result.mag.ln()).log10();
            if diff_to_next_slog < diff {
                result.mag = 1e10_f64.log10();
                result.layer += 1;
                diff -= diff_to_next_slog;
            }

            result.mag = result.mag.powf(10.0_f64.powf(diff));

            while subtract_layers_later > 0 {
                result.mag = result.mag.log10();
                subtract_layers_later -= 1;
            }
        }

        if diff < 0.0 {
            let mut subtract_layers_later: i64 = 0;

            while result.mag.is_finite() && result.mag < 10.0 {
                result.mag = 10.0_f64.powf(result.mag);
                subtract_layers_later += 1;
            }

            if result.mag > 1e10_f64 {
                result.mag = result.mag.log10();
                result.layer += 1;
            }

            let diff_to_next_slog = (1.0 / result.mag.log10()).log10();
            if diff_to_next_slog > diff {
                result.mag = 1e10_f64;
                result.layer -= 1;
                diff -= diff_to_next_slog;
            }

            result.mag = result.mag.powf(10.0_f64.powf(diff));

            while subtract_layers_later > 0 {
                result.mag = result.mag.log10();
                subtract_layers_later -= 1;
            }
        }

        while result.layer < 0 {
            result.layer += 1;
            result.mag = result.mag.log10();
        }

        result.normalize();
        result
    }

    /// Adds `diff` to the Decimal's slog(base) representation
    pub fn layer_add(&self, diff: f64, base: Decimal) -> Decimal {
        let slog_this = self.slog(Some(base)).to_number();
        let slog_dest = slog_this + diff;

        if slog_dest >= 0.0 {
            return base.tetrate(Some(slog_dest), None);
        }

        if !slog_dest.is_finite() {
            return Decimal::nan();
        }

        if slog_dest >= -1.0 {
            return base.tetrate(Some(slog_dest + 1.0), None).log(base);
        }

        base.tetrate(Some(slog_dest + 2.0), None)
            .log(base)
            .log(base)
    }

    /// Returns the product logarithm of the Decimal
    pub fn lambertw(&self) -> Result<Decimal, BreakEternityError> {
        if self < &Decimal::from_number(-0.3678794411710499) {
            return Err(BreakEternityError::LambertWError);
        }

        if self.mag < 0.0 {
            return Ok(Decimal::from_number(f_lambertw(self.to_number(), None)?));
        }

        if self.layer == 0 {
            return Ok(Decimal::from_number(f_lambertw(
                self.sign as f64 * self.mag,
                None,
            )?));
        }

        if self.layer == 1 || self.layer == 2 {
            return d_lambertw(*self, None);
        }

        Ok(Decimal::from_components_no_normalize(
            self.sign,
            self.layer - 1,
            self.mag,
        ))
    }

    /// Returns the super square root of the Decimal
    ///
    /// Essentially "what number, tetrated to height 2, equals this?"
    pub fn ssqrt(&self) -> Decimal {
        if self.sign == 1 && self.layer >= 3 {
            return Decimal::from_components_no_normalize(self.sign, self.layer - 1, self.mag);
        }

        let ln_x = self.ln();
        ln_x / ln_x.lambertw().expect("Expected number higher than -1")
    }

    /// The result of tetrating the Decimal `height` times in a row.
    pub fn pentate(&self, height: Option<Number>, payload: Option<Decimal>) -> Decimal {
        let mut height = height.unwrap_or(2.0_f64);
        let mut payload =
            payload.unwrap_or_else(|| Decimal::from_components_no_normalize(1, 0, 1.0));

        let old_height = height;
        height = height.trunc();
        let fract_height = old_height - height;

        if fract_height != 0.0 {
            if payload == Decimal::one() {
                height += 1.0;
                payload = Decimal::from_number(fract_height);
            } else if *self == Decimal::from_number(10.0) {
                payload = payload.layer_add_10(Decimal::from_number(fract_height));
            } else {
                payload = payload.layer_add(fract_height, *self);
            }
        }

        for i in 0..height as i64 {
            payload = self.tetrate(Some(payload.to_number()), None);
            if !payload.mag.is_finite() {
                return payload;
            }
            if i > 10 {
                return payload;
            }
        }

        payload
    }

    /// Returns the sin of the Decimal
    pub fn sin(&self) -> Decimal {
        if self.mag < 0.0 {
            return *self;
        }

        if self.layer == 0 {
            return Decimal::from_number((self.sign as f64 * self.mag).sin());
        }

        Decimal::from_components_no_normalize(0, 0, 0.0)
    }

    /// Returns the cos of the Decimal
    pub fn cos(&self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::one();
        }

        if self.layer == 0 {
            return Decimal::from_number((self.sign as f64 * self.mag).cos());
        }

        Decimal::from_components_no_normalize(0, 0, 0.0)
    }

    /// Returns the tan of the Decimal
    pub fn tan(&self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::from_number((self.sign as f64 * self.mag).tan());
        }

        if self.layer == 0 {
            return Decimal::from_number((self.sign as f64 * self.mag).tan());
        }

        Decimal::from_components_no_normalize(0, 0, 0.0)
    }

    /// Returns the asin of the Decimal
    pub fn asin(&self) -> Decimal {
        if self.mag < 0.0 {
            return *self;
        }

        if self.layer == 0 {
            return Decimal::from_number((self.sign as f64 * self.mag).asin());
        }

        Decimal::nan()
    }

    /// Returns the acos of the Decimal
    pub fn acos(&self) -> Decimal {
        if self.mag < 0.0 {
            return Decimal::from_number(self.to_number().acos());
        }

        if self.layer == 0 {
            return Decimal::from_number((self.sign as f64 * self.mag).acos());
        }

        Decimal::nan()
    }

    /// Returns the atan of the Decimal
    pub fn atan(&self) -> Decimal {
        if self.mag < 0.0 {
            return *self;
        }

        if self.layer == 0 {
            return Decimal::from_number((self.sign as f64 * self.mag).atan());
        }

        Decimal::from_number(std::f64::INFINITY.atan())
    }

    /// Returns the sinh of the Decimal
    pub fn sinh(&self) -> Decimal {
        (self.exp() - (-*self).exp()) / Decimal::from_number(2.0)
    }

    /// Returns the cosh of the Decimal
    pub fn cosh(&self) -> Decimal {
        (self.exp() + (-*self).exp()) / Decimal::from_number(2.0)
    }

    /// Returns the tanh of the Decimal
    pub fn tanh(&self) -> Decimal {
        self.sinh() / self.cosh()
    }

    /// Returns the asinh of the Decimal
    pub fn asinh(&self) -> Decimal {
        (*self + (self.sqr() + Decimal::from_number(1.0)).sqrt()).ln()
    }

    /// Returns the acosh of the Decimal
    pub fn acosh(&self) -> Decimal {
        (*self + (self.sqr() - Decimal::from_number(1.0)).sqrt()).ln()
    }

    /// Returns the atanh of the Decimal
    pub fn atanh(&self) -> Decimal {
        if self.abs() >= Decimal::from_number(1.0) {
            return Decimal::nan();
        }

        (*self + Decimal::from_number(1.0))
            / (Decimal::from_number(1.0) - *self).ln()
            / Decimal::from_number(2.0)
    }
}

impl PartialEq for Decimal {
    fn eq(&self, other: &Self) -> bool {
        // Special edge cases for NaN and infinities
        if self.mag.is_nan() && other.mag.is_nan() {
            return true;
        }

        if (self.mag.is_infinite() && self.mag.is_sign_positive())
            && (other.mag.is_infinite() && other.mag.is_sign_positive())
        {
            return true;
        }

        if (self.mag.is_infinite() && self.mag.is_sign_negative())
            && (other.mag.is_infinite() && other.mag.is_sign_negative())
        {
            return true;
        }

        self.sign == other.sign && self.layer == other.layer && (self.mag - other.mag).abs() < 1e-10
    }
}

impl Eq for Decimal {}

impl PartialOrd for Decimal {
    #[allow(clippy::comparison_chain)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.sign > other.sign {
            return Some(std::cmp::Ordering::Greater);
        }

        if self.sign < other.sign {
            return Some(std::cmp::Ordering::Less);
        }

        let cmp_abs = self.cmpabs(other) * self.sign;
        if cmp_abs > 0 {
            Some(std::cmp::Ordering::Greater)
        } else if cmp_abs < 0 {
            Some(std::cmp::Ordering::Less)
        } else {
            Some(std::cmp::Ordering::Equal)
        }
    }
}

impl Ord for Decimal {
    #[allow(clippy::comparison_chain)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.sign > other.sign {
            return std::cmp::Ordering::Greater;
        }

        if self.sign < other.sign {
            return std::cmp::Ordering::Less;
        }

        let cmp_abs = self.cmpabs(other) * self.sign;
        if cmp_abs > 0 {
            std::cmp::Ordering::Greater
        } else if cmp_abs < 0 {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

impl Add<Decimal> for Decimal {
    type Output = Decimal;

    fn add(self, rhs: Decimal) -> Self::Output {
        if !self.mag.is_finite() {
            return self;
        }

        if self.sign == 0 {
            return rhs;
        }
        if rhs.sign == 0 {
            return self;
        }

        if self.sign == -(rhs.sign) && self.layer == rhs.layer && (self.mag - rhs.mag).abs() < 1e-10
        {
            return Decimal::zero();
        }

        let a: Decimal;
        let b: Decimal;

        if self.layer >= 2 || rhs.layer >= 2 {
            return self.maxabs(rhs);
        }

        if self.cmpabs(&rhs) > 0 {
            a = self;
            b = rhs;
        } else {
            a = rhs;
            b = self;
        }

        if a.layer == 0 && b.layer == 0 {
            return Decimal::from_number(a.sign as f64 * a.mag + b.sign as f64 * b.mag);
        }

        let layer_a = a.layer * sign(a.mag) as i64;
        let layer_b = b.layer * sign(b.mag) as i64;

        if layer_a - layer_b >= 2 {
            return a;
        }

        if layer_a == 0 && layer_b == -1 {
            if (b.mag - a.mag.log10()).abs() > MAX_FLOAT_PRECISION as f64 {
                return a;
            } else {
                let mag_diff = 10.0_f64.powf(a.mag.log10() - b.mag);
                let mantissa = b.sign as f64 + (a.sign as f64 * mag_diff);
                return Decimal::from_components(sign(mantissa), 1, b.mag + mantissa.abs().log10());
            }
        }

        if layer_a == 1 && layer_b == 0 {
            if (a.mag - b.mag.log10()).abs() > MAX_FLOAT_PRECISION as f64 {
                return a;
            } else {
                let mag_diff = 10.0_f64.powf(a.mag - b.mag.log10());
                let mantissa = b.sign as f64 + (a.sign as f64 * mag_diff);
                return Decimal::from_components(
                    sign(mantissa),
                    1,
                    b.mag.log10() + mantissa.abs().log10(),
                );
            }
        }

        if (a.mag - b.mag).abs() > MAX_FLOAT_PRECISION as f64 {
            return a;
        }

        let mag_diff = 10.0_f64.powf(a.mag - b.mag);
        let mantissa = b.sign as f64 + (a.sign as f64 * mag_diff);
        let new_mag = b.mag + mantissa.abs().log10();
        Decimal::from_components(sign(mantissa), 1, new_mag)
    }
}

impl Sub<Decimal> for Decimal {
    type Output = Decimal;

    fn sub(self, rhs: Decimal) -> Self::Output {
        self + -rhs
    }
}

impl Mul<Decimal> for Decimal {
    type Output = Decimal;

    fn mul(self, rhs: Decimal) -> Self::Output {
        if self.sign == 0 || rhs.sign == 0 {
            return Decimal::zero();
        }

        if self.layer == rhs.layer && (self.mag - -rhs.mag).abs() < 1e-10 {
            return Decimal::from_components_no_normalize(self.sign * rhs.sign, 0, 1.0);
        }

        let a: Decimal;
        let b: Decimal;

        if (self.layer > rhs.layer) || (self.layer == rhs.layer && self.mag.abs() > rhs.mag.abs()) {
            a = self;
            b = rhs;
        } else {
            a = rhs;
            b = self;
        }

        if a.layer == 0 && b.layer == 0 {
            return Decimal::from_number(
                a.sign as f64 * b.sign as f64 * a.mag as f64 * b.mag as f64,
            );
        }

        if a.layer >= 3 || (a.layer - b.layer >= 2) {
            return Decimal::from_components(a.sign * b.sign, a.layer, a.mag);
        }

        if a.layer == 1 && b.layer == 0 {
            return Decimal::from_components(a.sign * b.sign, 1, a.mag + b.mag.log10());
        }

        if a.layer == 1 && b.layer == 1 {
            return Decimal::from_components(a.sign * b.sign, 1, a.mag + b.mag);
        }

        if a.layer == 2 && b.layer == 1 {
            let new_mag = Decimal::from_components(sign(a.mag), a.layer - 1, a.mag.abs())
                + Decimal::from_components(sign(b.mag), b.layer - 1, b.mag.abs());
            return Decimal::from_components(
                a.sign * b.sign,
                new_mag.layer + 1,
                new_mag.sign as f64 * new_mag.mag,
            );
        }

        if a.layer == 2 && b.layer == 2 {
            let new_mag = Decimal::from_components(sign(a.mag), a.layer - 1, a.mag.abs())
                + Decimal::from_components(sign(b.mag), b.layer - 1, b.mag.abs());
            return Decimal::from_components(
                a.sign * b.sign,
                new_mag.layer + 1,
                new_mag.sign as f64 * new_mag.mag,
            );
        }

        Decimal::inf()
    }
}

impl Div<Decimal> for Decimal {
    type Output = Decimal;

    /// Division of two decimals by multiplying the denominator by the reciprocal of the numerator.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Decimal) -> Self::Output {
        self * rhs.recip()
    }
}

impl AddAssign<Decimal> for Decimal {
    fn add_assign(&mut self, rhs: Decimal) {
        *self = *self + rhs;
    }
}

impl SubAssign<Decimal> for Decimal {
    fn sub_assign(&mut self, rhs: Decimal) {
        *self = *self - rhs;
    }
}

impl MulAssign<Decimal> for Decimal {
    fn mul_assign(&mut self, rhs: Decimal) {
        *self = *self * rhs;
    }
}

impl DivAssign<Decimal> for Decimal {
    fn div_assign(&mut self, rhs: Decimal) {
        *self = *self / rhs;
    }
}

impl RemAssign<Decimal> for Decimal {
    fn rem_assign(&mut self, rhs: Decimal) {
        *self = *self % rhs;
    }
}

impl Rem<Decimal> for Decimal {
    type Output = Decimal;

    fn rem(self, rhs: Decimal) -> Self::Output {
        if rhs == Decimal::zero() {
            return Decimal::zero();
        }

        if self.sign * rhs.sign == -1 {
            return self.abs().rem(rhs.abs()).neg();
        }

        if self.sign == -1 {
            return self.abs().rem(rhs.abs());
        }

        self - (self / rhs).floor() * rhs
    }
}

impl Neg for Decimal {
    type Output = Decimal;

    fn neg(self) -> Decimal {
        Decimal::from_components_no_normalize(-self.sign, self.layer, self.mag)
    }
}

#[cfg(feature = "serde")]
mod serde {
    use super::*;

    impl serde_crate::Serialize for Decimal {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: serde_crate::Serializer,
        {
            serializer.serialize_str(&self.to_string())
        }
    }

    impl<'de> serde_crate::Deserialize<'de> for Decimal {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: serde_crate::Deserializer<'de>,
        {
            let d = String::deserialize(deserializer)?;
            let dec: Result<Decimal, BreakEternityError> = d.as_str().try_into();
            dec.map_err(|_| serde_crate::de::Error::custom("Could not parse Decimal"))
        }
    }
}


impl LowerExp for Decimal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if *self == Decimal::inf() {
            return write!(f, "Infinity");
        }

        if *self == Decimal::neg_inf() {
            return write!(f, "-Infinity");
        }

        if *self == Decimal::nan() {
            return write!(f, "NaN");
        }

        let precision = f.precision().unwrap_or(2);

        if self.layer == 0 {
            return write!(f, "{:.*e}", precision, self.sign as f64 * self.mag);
        }
        write!(
            f,
            "{}",
            self.to_string_with_decimal_places(precision, Some(true))
        )
    }
}

impl UpperExp for Decimal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if *self == Decimal::inf() {
            return write!(f, "Infinity");
        }

        if *self == Decimal::neg_inf() {
            return write!(f, "-Infinity");
        }

        if *self == Decimal::nan() {
            return write!(f, "NaN");
        }

        let precision = f.precision().unwrap_or(2);

        if self.layer == 0 {
            return write!(f, "{:.*E}", precision, self.sign as f64 * self.mag);
        }
        write!(
            f,
            "{}",
            self.to_string_with_decimal_places(precision, Some(false))
        )
    }
}

impl Display for Decimal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if *self == Decimal::inf() {
            return write!(f, "Infinity");
        }

        if *self == Decimal::neg_inf() {
            return write!(f, "-Infinity");
        }

        if *self == Decimal::nan() {
            return write!(f, "NaN");
        }

        if self.layer == 0 {
            if (self.mag < 1e21 && self.mag > 1e-7) || self.mag == 0.0 {
                return write!(f, "{}", self.sign as f64 * self.mag);
            }
            return write!(f, "{}e{}", self.m(), self.e());
        }
        if self.layer == 1 {
            return write!(f, "{}e{}", self.m(), self.e());
        }
        if self.layer <= MAX_ES_IN_A_ROW as i64 {
            return write!(
                f,
                "{}{}{}",
                if self.sign == -1 { "-" } else { "" },
                "e".repeat(self.layer as usize),
                self.mag
            );
        }
        write!(
            f,
            "{}(e^{}){}",
            if self.sign == -1 { "-" } else { "" },
            self.layer,
            self.mag
        )
    }
}

macro_rules! impl_from_primitive {
    ($prim_type:ty) => {
        impl From<$prim_type> for Decimal {
            fn from(prim: $prim_type) -> Self {
                Decimal::from_number(prim as f64)
            }
        }
    };
}

macro_rules! impl_ops_primitive {
    ($prim_type:ty) => {
        impl Add<$prim_type> for Decimal {
            type Output = Decimal;

            fn add(self, rhs: $prim_type) -> Self::Output {
                self + Decimal::from_number(rhs as f64)
            }
        }

        impl Sub<$prim_type> for Decimal {
            type Output = Decimal;

            fn sub(self, rhs: $prim_type) -> Self::Output {
                self - Decimal::from_number(rhs as f64)
            }
        }

        impl Mul<$prim_type> for Decimal {
            type Output = Decimal;

            fn mul(self, rhs: $prim_type) -> Self::Output {
                self * Decimal::from_number(rhs as f64)
            }
        }

        impl Div<$prim_type> for Decimal {
            type Output = Decimal;

            fn div(self, rhs: $prim_type) -> Self::Output {
                self / Decimal::from_number(rhs as f64)
            }
        }

        impl Rem<$prim_type> for Decimal {
            type Output = Decimal;

            fn rem(self, rhs: $prim_type) -> Self::Output {
                self % Decimal::from_number(rhs as f64)
            }
        }

        impl AddAssign<$prim_type> for Decimal {
            fn add_assign(&mut self, rhs: $prim_type) {
                *self = *self + rhs;
            }
        }

        impl SubAssign<$prim_type> for Decimal {
            fn sub_assign(&mut self, rhs: $prim_type) {
                *self = *self - rhs;
            }
        }

        impl MulAssign<$prim_type> for Decimal {
            fn mul_assign(&mut self, rhs: $prim_type) {
                *self = *self * rhs;
            }
        }

        impl DivAssign<$prim_type> for Decimal {
            fn div_assign(&mut self, rhs: $prim_type) {
                *self = *self / rhs;
            }
        }

        impl RemAssign<$prim_type> for Decimal {
            fn rem_assign(&mut self, rhs: $prim_type) {
                *self = *self % rhs;
            }
        }
    };
}

impl_from_primitive!(i8);
impl_from_primitive!(i16);
impl_from_primitive!(i32);
impl_from_primitive!(i64);
impl_from_primitive!(u8);
impl_from_primitive!(u16);
impl_from_primitive!(u32);
impl_from_primitive!(u64);
impl_from_primitive!(f32);
impl_from_primitive!(f64);

impl_ops_primitive!(i8);
impl_ops_primitive!(i16);
impl_ops_primitive!(i32);
impl_ops_primitive!(i64);
impl_ops_primitive!(u8);
impl_ops_primitive!(u16);
impl_ops_primitive!(u32);
impl_ops_primitive!(u64);
impl_ops_primitive!(f32);
impl_ops_primitive!(f64);

impl TryFrom<&str> for Decimal {
    type Error = BreakEternityError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut value = s.to_string();
        if *IGNORE_COMMAS {
            value = value.replace(",", "");
        } else if *COMMAS_ARE_DECIMAL_POINTS {
            value = value.replace(",", ".");
        }
        let value = value.as_str();

        let pentation_parts: Vec<&str> = value.split("^^^").collect();
        if pentation_parts.len() == 2 {
            let base = pentation_parts[0].parse::<Number>();
            if let Err(parse_error) = base {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let base = base.unwrap();
            let height = pentation_parts[1].parse::<Number>();
            if let Err(parse_error) = height {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let height = height.unwrap();
            let mut payload = 1.0;
            let height_parts = pentation_parts[1].split(';').collect::<Vec<&str>>();
            if height_parts.len() == 2 {
                let payload_parsed = height_parts[1].parse::<Number>();
                if let Err(parse_error) = payload_parsed {
                    return Err(BreakEternityError::ParseError {
                        parsed: s.to_string(),
                        error: parse_error,
                    });
                }
                payload = payload_parsed.unwrap();
                if !payload.is_finite() {
                    payload = 1.0;
                }
            }

            if base.is_finite() && height.is_finite() {
                return Ok(Decimal::from_number(base)
                    .pentate(Some(height), Some(Decimal::from_number(payload))));
            }
        }

        let tetration_parts: Vec<&str> = value.split("^^").collect();
        if tetration_parts.len() == 2 {
            let base = tetration_parts[0].parse::<Number>();
            if let Err(parse_error) = base {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let base = base.unwrap();
            let height = tetration_parts[1].parse::<Number>();
            if let Err(parse_error) = height {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let height = height.unwrap();
            let mut payload = 1.0;
            let height_parts = tetration_parts[1].split(';').collect::<Vec<&str>>();
            if height_parts.len() == 2 {
                let payload_parsed = height_parts[1].parse::<Number>();
                if let Err(parse_error) = payload_parsed {
                    return Err(BreakEternityError::ParseError {
                        parsed: s.to_string(),
                        error: parse_error,
                    });
                }
                payload = payload_parsed.unwrap();
                if !payload.is_finite() {
                    payload = 1.0;
                }
            }

            if base.is_finite() && height.is_finite() {
                return Ok(Decimal::from_number(base)
                    .tetrate(Some(height), Some(Decimal::from_number(payload))));
            }
        }

        let pow_parts = value.split('^').collect::<Vec<&str>>();
        if pow_parts.len() == 2 {
            let base = pow_parts[0].parse::<Number>();
            if let Err(parse_error) = base {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let base = base.unwrap();
            let exponent = pow_parts[1].parse::<Number>();
            if let Err(parse_error) = exponent {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let exponent = exponent.unwrap();
            if base.is_finite() && exponent.is_finite() {
                return Ok(Decimal::from_number(base).pow(Decimal::from_number(exponent)));
            }
        }

        let value = value.trim().to_lowercase();
        let value = value.as_str();

        let mut pt_parts = value.split("pt").collect::<Vec<&str>>();
        if pt_parts.len() == 2 {
            let base: f64 = 10.0;
            let height = pt_parts[0].parse::<Number>();
            if let Err(parse_error) = height {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let height = height.unwrap();
            let tmp = pt_parts[1].replace("(", "").replace(")", "");
            pt_parts[1] = tmp.as_str();

            let payload = pt_parts[1].parse::<Number>();
            if let Err(parse_error) = payload {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let mut payload = payload.unwrap();
            if !payload.is_finite() {
                payload = 1.0;
            }
            if base.is_finite() && height.is_finite() {
                // tetrate again
                return Ok(Decimal::from_number(base)
                    .tetrate(Some(height), Some(Decimal::from_number(payload))));
            }
        }

        let mut p_parts = value.split('p').collect::<Vec<&str>>();
        if p_parts.len() == 2 {
            let base: f64 = 10.0;
            let height = p_parts[0].parse::<Number>();
            if let Err(parse_error) = height {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let height = height.unwrap();
            let tmp = p_parts[1].replace("(", "").replace(")", "");
            p_parts[1] = tmp.as_str();

            let payload = p_parts[1].parse::<Number>();
            if let Err(parse_error) = payload {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let mut payload = payload.unwrap();
            if !payload.is_finite() {
                payload = 1.0;
            }
            if base.is_finite() && height.is_finite() {
                // another tetrate
                return Ok(Decimal::from_number(base)
                    .tetrate(Some(height), Some(Decimal::from_number(payload))));
            }
        }

        let e_parts = value.split('e').collect::<Vec<&str>>();
        let e_count = e_parts.len() - 1;

        if e_count == 0 {
            let number_attempt = value.parse::<Number>();
            if let Err(parse_error) = number_attempt {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let number_attempt = number_attempt.unwrap();
            if number_attempt.is_finite() {
                return Ok(Decimal::default().set_from_number(number_attempt));
            }
        } else if e_count == 1 {
            let number_attempt = value.parse::<Number>();
            if let Err(parse_error) = number_attempt {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let number_attempt = number_attempt.unwrap();
            if number_attempt.is_finite() && number_attempt != 0.0 {
                return Ok(Decimal::default().set_from_number(number_attempt));
            }
        }

        let new_parts = value.split("e^").collect::<Vec<&str>>();
        if new_parts.len() == 2 {
            let mut dec = Decimal {
                sign: 1,
                ..Default::default()
            };
            if new_parts[0].starts_with('-') {
                dec.sign = -1;
            }

            let mut layer_string = "".to_string();
            for (i, c) in new_parts[1].chars().enumerate() {
                if c.is_numeric()
                    || c == '+'
                    || c == '-'
                    || c == '.'
                    || c == 'e'
                    || c == ','
                    || c == '/'
                {
                    layer_string.push(c);
                } else {
                    let layer = layer_string.parse::<Number>();
                    if let Err(parse_error) = layer {
                        return Err(BreakEternityError::ParseError {
                            parsed: s.to_string(),
                            error: parse_error,
                        });
                    }
                    dec.layer = layer.unwrap() as i64;
                    let mag = new_parts[1][i + 1..].parse::<Number>();
                    if let Err(parse_error) = mag {
                        return Err(BreakEternityError::ParseError {
                            parsed: s.to_string(),
                            error: parse_error,
                        });
                    }
                    let mag = mag.unwrap();
                    dec.mag = mag;
                    dec.normalize();
                    return Ok(dec);
                }
            }
        }

        let mut dec = Decimal::default();

        if e_count < 1 {
            return Ok(dec);
        }
        let mantissa = e_parts[0].parse::<Number>();
        if let Err(parse_error) = mantissa {
            return Err(BreakEternityError::ParseError {
                parsed: s.to_string(),
                error: parse_error,
            });
        }
        let mantissa = mantissa.unwrap();
        if mantissa == 0.0 {
            return Ok(dec);
        }
        let exponent = e_parts.last().unwrap().parse::<Number>();
        if let Err(parse_error) = exponent {
            return Err(BreakEternityError::ParseError {
                parsed: s.to_string(),
                error: parse_error,
            });
        }
        let mut exponent = exponent.unwrap();
        if e_count >= 2 {
            let me = e_parts[e_parts.len() - 2].parse::<Number>();
            if let Err(parse_error) = me {
                return Err(BreakEternityError::ParseError {
                    parsed: s.to_string(),
                    error: parse_error,
                });
            }
            let me = me.unwrap();
            if me.is_finite() {
                exponent *= sign(me) as f64;
                exponent += f_maglog10(me);
            }
        }

        if !mantissa.is_finite() {
            dec.sign = if e_parts[0] == "-" { -1 } else { 1 };
            dec.layer = e_count as i64;
            dec.mag = exponent;
        } else if e_count == 1 {
            dec.sign = sign(mantissa);
            dec.layer = 1;
            dec.mag = exponent + mantissa.abs().log10();
        } else {
            dec.sign = sign(mantissa);
            dec.layer = e_count as i64;
            if e_count == 2 {
                return Ok(
                    Decimal::from_components(1, 2, exponent) * Decimal::from_number(mantissa)
                );
            } else {
                // mantissa is way too small
                dec.mag = exponent;
            }
        }

        dec.normalize();
        Ok(dec)
    }
}
