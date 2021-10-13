use super::Decimal;

#[test]
fn decimal() {
    assert_eq!(Decimal::from_number(0.0).to_string(), "0");
	assert_eq!(Decimal::from_number(f64::NAN).to_string(), "NaN");
	assert_eq!(Decimal::from_number(f64::INFINITY).to_string(), "Infinity");
	assert_eq!(Decimal::from_number(f64::NEG_INFINITY).to_string(), "-Infinity");

	assert_eq!(Decimal::from_number(100.0).to_string(), "100");
	assert_eq!(Decimal::from_number(1e12).to_string(), "1000000000000");
	assert_eq!(Decimal::from_number(1.79e3).to_string(), "1790");
	assert_eq!(Decimal::from_number(1e308).to_string(), "1e308");
}

#[test]
fn simple_maths() {
	let a = Decimal::from_number(4.0);
	let b = Decimal::from_number(2.0);

	assert_eq!(a + b, Decimal::from_number(6.0));
	assert_eq!(a - b, Decimal::from_number(2.0));
	assert_eq!(a * b, Decimal::from_number(8.0));
	assert_eq!(a / b, Decimal::from_number(2.0));
}

#[test]
fn ops() {
    let a = Decimal::from_mantissa_exponent_no_normalize(3.224, 54.0);
    let b = Decimal::from_mantissa_exponent_no_normalize(1.24, 53.0);
    let c = Decimal::from_mantissa_exponent_no_normalize(3.1, 52.0);

    assert_eq!(a + b, Decimal::from_mantissa_exponent_no_normalize(3.348, 54.0));
	assert_eq!(a - b, Decimal::from_mantissa_exponent_no_normalize(3.1, 54.0));
	assert_eq!(a * b, Decimal::from_mantissa_exponent_no_normalize(3.9977600000000004, 107.0));
	assert_eq!(a / b, Decimal::from_mantissa_exponent_no_normalize(2.6, 1.0));

	assert_eq!(a + c, Decimal::from_mantissa_exponent_no_normalize(3.255, 54.0));
	assert_eq!(a - c, Decimal::from_mantissa_exponent_no_normalize(3.193, 54.0));
	assert_eq!(a * c, Decimal::from_mantissa_exponent_no_normalize(9.9944, 106.0));
	assert_eq!(a / c, Decimal::from_mantissa_exponent_no_normalize(1.04, 2.0));

	assert_eq!(b + c, Decimal::from_mantissa_exponent_no_normalize(1.55, 53.0));
	assert_eq!(b - c, Decimal::from_mantissa_exponent_no_normalize(9.3, 52.0));
	assert_eq!(b * c, Decimal::from_mantissa_exponent_no_normalize(3.844, 105.0));
	assert_eq!(b / c, Decimal::from_mantissa_exponent_no_normalize(3.9999999999999996, 0.0));
}

#[test]
fn rem() {
    let a = Decimal::from_number(5.0);
    let b = Decimal::from_number(2.0);

    assert_eq!(a % b, Decimal::from_number(1.0));
}

#[test]
#[allow(clippy::bool_assert_comparison)]
fn cmp() {
	let a = Decimal::from_mantissa_exponent_no_normalize(3.224, 54.0);
	let b = Decimal::from_mantissa_exponent_no_normalize(1.24, 53.0);
	let c = Decimal::from_mantissa_exponent_no_normalize(3.1, 52.0);
	let d = Decimal::from_mantissa_exponent_no_normalize(3.224, 54.0);

	assert_eq!(a == b, false);
	assert_eq!(a == d, true);
	assert_eq!(b == d, false);

	assert_eq!(a >= b, true);
	assert_eq!(a >= d, true);
	assert_eq!(b >= d, false);

	assert_eq!(a > b, true);
	assert_eq!(a > d, false);
	assert_eq!(b > d, false);

	assert_eq!(a <= b, false);
	assert_eq!(a <= d, true);
	assert_eq!(b <= d, true);

	assert_eq!(a < b, false);
	assert_eq!(a < d, false);
	assert_eq!(b < d, true);

	assert_eq!(a.max(b), a);
	assert_eq!(a.max(c), a);
	assert_eq!(b.max(c), b);

	assert_eq!(a.min(b), b);
	assert_eq!(a.min(c), c);
	assert_eq!(b.min(c), c);

	assert_eq!(a.clamp(c, b), b);
	assert_eq!(b.clamp(c, a), b);
	assert_eq!(c.clamp(b, b), b);
}

#[test]
fn neg_abs() {
    assert_eq!(-Decimal::from_number(456.7), Decimal::from_mantissa_exponent_no_normalize(-4.567, 2.0));
	assert_eq!(-Decimal::from_number(1.23e48), Decimal::from_mantissa_exponent_no_normalize(-1.23, 48.0));

	assert_eq!(Decimal::from_number(-456.7).abs(), Decimal::from_mantissa_exponent_no_normalize(4.567, 2.0));
	assert_eq!(Decimal::from_number(-1.23e48).abs(), Decimal::from_mantissa_exponent_no_normalize(1.23, 48.0));
}