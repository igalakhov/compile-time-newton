#include <iostream>

using ll = long long;

template <ll, ll = 1> struct Rational;

template <typename, typename> struct Complex;

template <typename...> struct Polynomial;

// simplification of rationals

template <typename, ll = 1000000> struct is_thresh;
template <ll n, ll d, ll L> struct is_thresh<Rational<n, d>, L> {
  static constexpr bool value =
      ((n > 0 ? n : -n) <= L) && ((d > 0 ? d : -d) <= L);
};

template <typename, ll> struct threshhold;
template <typename R> using threshold_t = typename threshhold<R, 1000000>::type;

template <ll n, ll d, ll L> struct threshhold<Rational<n, d>, L> {
  using type =
      typename std::conditional_t<is_thresh<Rational<n, d>, L>::value,
                                  std::type_identity<Rational<n, d>>,
                                  threshhold<Rational<n / 2, d / 2>, L>>::type;
};

template <ll, ll, typename = void> struct Gcd;

template <ll a, ll b> struct Gcd<a, b, std::enable_if_t<a >= 0 && (b > 0)>> {
  static constexpr auto value = Gcd<b, a % b>::value;
};

template <ll a, ll b> struct Gcd < a, b, std::enable_if_t<a<0 && (b > 0)>> {
  static constexpr auto value = -Gcd<-a, b>::value;
};

template <ll a, ll b> struct Gcd<a, b, std::enable_if_t<a >= 0 && (b < 0)>> {
  static constexpr auto value = -Gcd<a, -b>::value;
};

template <ll a, ll b> struct Gcd < a, b, std::enable_if_t<a<0 && (b < 0)>> {
  static constexpr auto value = Gcd<a, -b>::value;
};

template <ll a> struct Gcd<a, 0> {
  static constexpr auto value = a;
};

template <typename> struct Simplify;
template <typename R> using simplify_t = typename Simplify<R>::type;

template <ll n, ll d> struct Simplify<Rational<n, d>> {
  using type = Rational<n / Gcd<n, d>::value, d / Gcd<n, d>::value>;
};

// conversion to complex

template <typename T> struct to_complex;
template <typename T> using to_complex_t = typename to_complex<T>::type;

template <ll n, ll d> struct to_complex<Rational<n, d>> {
  using type = Complex<Rational<n, d>, Rational<0>>;
};

template <typename R, typename I> struct to_complex<Complex<R, I>> {
  using type = Complex<R, I>;
};

// addition

template <typename, typename> struct stupid_add;
template <ll n1, ll d1, ll n2, ll d2>
struct stupid_add<Rational<n1, d1>, Rational<n2, d2>> {
  using type = simplify_t<Rational<n1 * d2 + d1 * n2, d1 * d2>>;
};

template <typename, typename> struct Add;
template <typename V1, typename V2> using add_t = typename Add<V1, V2>::type;

template <ll n1, ll d1, ll n2, ll d2>
struct Add<Rational<n1, d1>, Rational<n2, d2>> {
  using type = typename std::conditional_t<
      is_thresh<Rational<n1, d1>>::value && is_thresh<Rational<n2, d2>>::value,
      stupid_add<Rational<n1, d1>, Rational<n2, d2>>,
      Add<threshold_t<Rational<n1, d1>>, threshold_t<Rational<n2, d2>>>>::type;
};

template <typename R1, typename I1, typename R2, typename I2>
struct Add<Complex<R1, I1>, Complex<R2, I2>> {
  using type = Complex<add_t<R1, R2>, add_t<I1, I2>>;
};

template <typename V1, typename V2> struct Add {
  using type = add_t<to_complex_t<V1>, to_complex_t<V2>>;
};

// multiplication

template <typename, typename> struct stupid_mult;
template <ll n1, ll d1, ll n2, ll d2>
struct stupid_mult<Rational<n1, d1>, Rational<n2, d2>> {
  using type = simplify_t<Rational<n1 * n2, d1 * d2>>;
};

template <typename, typename> struct Mult;
template <typename V1, typename V2> using mult_t = typename Mult<V1, V2>::type;

template <ll n1, ll d1, ll n2, ll d2>
struct Mult<Rational<n1, d1>, Rational<n2, d2>> {
  using type = typename std::conditional_t<
      is_thresh<Rational<n1, d1>>::value && is_thresh<Rational<n2, d2>>::value,
      stupid_mult<Rational<n1, d1>, Rational<n2, d2>>,
      Mult<threshold_t<Rational<n1, d1>>, threshold_t<Rational<n2, d2>>>>::type;
};

template <typename R1, typename I1, typename R2, typename I2>
struct Mult<Complex<R1, I1>, Complex<R2, I2>> {
  using type =
      Complex<add_t<mult_t<R1, R2>, mult_t<Rational<-1>, mult_t<I1, I2>>>,
              add_t<mult_t<R1, I2>, mult_t<R2, I1>>>;
};

template <typename V1, typename V2> struct Mult {
  using type = mult_t<to_complex_t<V1>, to_complex_t<V2>>;
};

// division
template <typename, typename> struct Divide;
template <typename V1, typename V2>
using divide_t = typename Divide<V1, V2>::type;

template <ll n1, ll d1, ll n2, ll d2>
struct Divide<Rational<n1, d1>, Rational<n2, d2>> {
  using type = mult_t<Rational<n1, d1>, Rational<d2, n2>>;
};

template <ll n, ll d, typename R, typename I>
struct Divide<Complex<R, I>, Rational<n, d>> {
  using type =
      Complex<divide_t<R, Rational<n, d>>, divide_t<I, Rational<n, d>>>;
};

template <typename R1, typename I1, typename R2, typename I2>
struct Divide<Complex<R1, I1>, Complex<R2, I2>> {
  using type =
      divide_t<mult_t<Complex<R1, I1>, Complex<R2, mult_t<Rational<-1>, I2>>>,
               add_t<mult_t<R2, R2>, mult_t<I2, I2>>>;
};

template <typename V1, typename V2> struct Divide {
  using type = divide_t<to_complex_t<V1>, to_complex_t<V2>>;
};

// polynomial operations

template <typename> struct degree;

template <typename... Cp> struct degree<Polynomial<Cp...>> {
  static constexpr ll value = sizeof...(Cp);
};

template <typename, typename> struct evaluate_polynomial;

template <typename P, typename V>
using evaluate_polynomial_t = typename evaluate_polynomial<P, V>::type;

template <typename V> struct evaluate_polynomial<Polynomial<>, V> {
  using type = Rational<0>;
};

template <typename C0, typename... Cr, typename V>
struct evaluate_polynomial<Polynomial<C0, Cr...>, V> {
  using type =
      add_t<mult_t<evaluate_polynomial_t<Polynomial<Cr...>, V>, V>, C0>;
};

template <typename, ll> struct derivative;
template <typename P>
using derivative_t = typename derivative<P, degree<P>::value>::type;

template <ll N> struct derivative<Polynomial<>, N> {
  using type = Polynomial<>;
};

template <typename C0, typename... Cp, ll N>
struct derivative<Polynomial<C0, Cp...>, N> {

  template <typename, typename...> struct smash;
  template <typename F, typename... R>
  using smash_t = typename smash<F, R...>::type;

  template <typename F, typename... R> struct smash<F, Polynomial<R...>> {
    using type = Polynomial<F, R...>;
  };

  using type = typename std::conditional_t<
      degree<Polynomial<C0, Cp...>>::value == N,
      derivative<Polynomial<Cp...>, N>,
      smash<mult_t<C0, Rational<N - degree<Polynomial<Cp...>>::value - 1>>,
            typename derivative<Polynomial<Cp...>, N>::type>>::type;
};

// newton map
template <typename> struct NewtonMap;

template <typename P> struct NewtonMap {
  using Pd = derivative_t<P>;

  template <typename V> struct evaluate {
    using type =
        add_t<V, mult_t<Rational<-1>, divide_t<evaluate_polynomial_t<P, V>,
                                               evaluate_polynomial_t<Pd, V>>>>;
  };

  template <typename V> using evaluate_t = typename evaluate<V>::type;
};

// applying many times
template <typename, size_t, typename> struct iterate;
template <typename Map, size_t N, typename V>
using iterate_t = typename iterate<Map, N, V>::type;

template <typename Map, size_t N, typename V> struct iterate {
  using type = typename std::conditional_t<
      N == 0, std::type_identity<V>,
      iterate<Map, N - 1, typename Map::template evaluate_t<V>>>::type;
};

// utils

template <typename> struct Printer;
template <ll n, ll d> struct Printer<Rational<n, d>> {
  static void print(const char *end = "\n") {
    std::cout << n << "/" << d << end;
  }
};

template <typename R, typename I> struct Printer<Complex<R, I>> {
  static void print(const char *end = "\n") {
    Printer<R>::print(" + ");
    Printer<I>::print("i\n");
  }
};

int main() {
  using P = Polynomial<Rational<1>, Rational<0>, Rational<1>>;
  using Map = NewtonMap<P>;
  using V = iterate_t<Map, 100, Complex<Rational<1>, Rational<1>>>;

  Printer<V>::print();
}