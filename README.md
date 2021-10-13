# break-eternity

A port of [Patashu's break_eternity.js](https://github.com/Patashu/break_eternity.js).

A numerical library to represent numbers as large as 10^^1e308 and as 'small' as 10^-(10^^1e308).

This library focuses less on precision and more on speed. It is intended to be used by games.

## Additional Features

This crate has 2 more features that can be used:

* `serde`, which adds support for serde
* `godot`, which adds `FromVariant` and `ToVariant` from `gdnative` to the derived traits

By default, both features are disabled. If you want this library to include support for a different library, please open an issue and tell me about it. I would be more than happy to add more support for game engines, since this is a library that's supposed to be used for games.

## Internal Representation

The internal representation of a Decimal is as follows:
`Decimal::from_components(sign, layer, mag)` == `sign * 10^10^10^...(layer times) mag`. So a layer 0 number is just `sign * mag`, a layer 1 number is `sign * 10^mag`, a layer 2 number is `sign * 10^10^mag` and so on.

If `layer > 0` and `mag < 0.0`, then the number's exponent is negative, e.g. `sign * 10^-10^10^10^ ... mag`.

* `sign` is -1, 0 or 1
* `layer` is a non-negative integer
* `mag` is an `f64`, normalized as follows: If it is above 1e15, `log10(mag)` it and increment layer. If it is below `log10(9e15)` (about 15.954) and `layer > 0`, `10.0_f64.powf(mag)` it and decrement layer. At layer 0, sign is extracted from negative mags. Zeroes (`sign == 0 || (mag == 0.0 && layer == 0)`) become `0, 0, 0` in all fields.

Decimal implements `Copy` and `Clone`, so it can be safely dereferenced without a fuss.

## Creating a Decimal

You can create a Decimal using `Decimal::from_number(f64)`, `Decimal::try_from(&str)` or manually using `Decimal::from_components(sign, layer, mag)` or `Decimal::from_mantissa_exponent(mantissa, exponent)`.

If you use the struct initialization syntax, please make sure to run the `normalize()` function to normalize the Decimal.

### Accepted String representations

```plain
M === M
eX === 10^X
MeX === M*10^X
eXeY === 10^(XeY)
MeXeY === M*10^(XeY)
eeX === 10^10^X
eeXeY === 10^10^(XeY)
eeeX === 10^10^10^X
eeeXeY === 10^10^10^(XeY)
eeee... (N es) X === 10^10^10^ ... (N 10^s) X
(e^N)X === 10^10^10^ ... (N 10^s) X
N PT X === 10^10^10^ ... (N 10^s) X
N PT (X) === 10^10^10^ ... (N 10^s) X
NpX === 10^10^10^ ... (N 10^s) X
X^Y === X^Y
X^^N === X^X^X^ ... (N X^s) 1
X^^N;Y === X^X^X^ ... (N X^s) Y
X^^^N === X^^X^^X^^ ... (N X^^s) 1
X^^^N;Y === X^^X^^X^^ ... (N X^^s) Y
```

## Operations

Thanks to the power of Rust traits, you can simply use the regular operators (`+`, `-`, `*`, `/`, `%`, `+=`, `-=`, `*=`, `/=`, `%=`) for math operations as well as other mathematical functions, such as: `abs, neg, round, floor, ceil, trunc, recip, cmp, cmpabs, max, min, maxabs, minabs, log, log10, ln, pow, root, factorial, gamma, exp, sqrt, tetrate, iteratedexp, iteratedlog, layer_add_10, layer_add, slog, ssqrt, lambertw, pentate`.

Equality is handled in a special way, such that if both sides are NaN or Infinity, they are equal. Other than that, Decimals are considered equal to a precision of `1e-10`.

As seen above, the modulo operator is also implemented properly and should be just as accurate as other operations.

Another simplicity feature is the implementation for said operators for all primitive number types, which means you can add like this `Decimal::from_number(1.0) + 2.0`. Conversions also work the same way, any primitve number type can be converted to Decimal using `from()` and `into()`.

## Note to bugs

Even though I ported this library, I am not very knowledgeable about math. In fact, I have no idea how to properly apply stuff like tetration, super log etc. So if there are any issues with these functions, you will probably have to explain me what I have to change  in the code.

I know this isn't very professional, but I lack the time to get a degree in math or further educate myself about the topic. This crate was pretty much born out of necessity, since I want to implement this in games I develop.

Other than that, all other kinds of bugs are appreciated.

## Afterword

This crate took a long time to port. Mostly because it IS a huge library. And being a single person, future ports of changes to `break_eternity.js` may also take a while.

Please consider contributing and actively opening pull requests, if you have improvements. It would help out a lot.

There definetly will be bugs, though I am more than willing to fix them.

Also, if you like this library, why not leave a star on this and [Patashu's break_eternity.js](https://github.com/Patashu/break_eternity.js)? It helps boost the popularity of our packages.

## Special Thanks

* Patashu (for taking the time to write this HUGE library)
* Naruyoko (for OmegaNum.js's modulo implementation, which I shamelessly copied)
