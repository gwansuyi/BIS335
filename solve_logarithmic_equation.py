"""
Find the smallest integer m (>= 2) from which 2*[ln(m)]^2 < m ALWAYS holds.

Equivalently, using the rearranged form (dividing by 2*m*ln(m), valid for m > 1):
    1/(2*ln(m)) > ln(m)/m

Both forms are compared side-by-side to verify they give the same answer.
"""

import math


def f(m):
    """f(m) = 2*[ln(m)]^2 - m"""
    return 2 * math.log(m) ** 2 - m


def f_prime(m):
    """f'(m) = 4*ln(m)/m - 1"""
    return 4 * math.log(m) / m - 1


def bisection(a, b, tol=1e-15, max_iter=200):
    """Find root of f in [a, b] by bisection."""
    for _ in range(max_iter):
        mid = (a + b) / 2
        if f(mid) == 0 or (b - a) / 2 < tol:
            return mid
        if f(a) * f(mid) < 0:
            b = mid
        else:
            a = mid
    return (a + b) / 2


def newton(m0, tol=1e-15, max_iter=100):
    """Find root of f using Newton's method starting from m0."""
    m = m0
    for _ in range(max_iter):
        fm = f(m)
        fpm = f_prime(m)
        if abs(fpm) < 1e-30:
            break
        m_new = m - fm / fpm
        if m_new <= 0:
            m_new = m / 2  # safeguard: stay positive
        if abs(m_new - m) < tol:
            return m_new
        m = m_new
    return m


if __name__ == "__main__":
    print("=" * 60)
    print("Compare two equivalent forms (for m >= 2):")
    print("  Form A:  2*[ln(m)]^2 < m")
    print("  Form B:  1/(2*ln(m)) > ln(m)/m")
    print("=" * 60)

    # --- Find the last root (the one that matters for "always holds") ---
    root3 = newton(bisection(13.0, 14.0))
    print(f"\nLast root of 2*[ln(m)]^2 = m:  m ≈ {root3:.10f}")
    print(f"=> For all m > {root3:.4f}, both inequalities hold permanently.")

    # --- Side-by-side comparison ---
    print(f"\n  {'m':>4s}   {'2*ln²(m)':>10s} {'<':>2s} {'m':>4s}?   "
          f"{'1/(2ln m)':>10s} {'> ':>2s} {'ln(m)/m':>10s}?   {'Same?':>5s}")
    print(f"  {'----':>4s}   {'----------':>10s} {' ':>2s} {'----':>4s}    "
          f"{'----------':>10s} {'  ':>2s} {'----------':>10s}    {'-----':>5s}")

    all_same = True
    first_always = None
    for m in range(2, 21):
        ln_m = math.log(m)

        # Form A: 2*[ln(m)]^2 < m
        lhs_a = 2 * ln_m ** 2
        form_a = lhs_a < m

        # Form B: 1/(2*ln(m)) > ln(m)/m
        lhs_b = 1.0 / (2 * ln_m)
        rhs_b = ln_m / m
        form_b = lhs_b > rhs_b

        same = form_a == form_b
        if not same:
            all_same = False

        a_str = "YES" if form_a else "no"
        b_str = "YES" if form_b else "no"
        s_str = "OK" if same else "DIFF"

        print(f"  {m:4d}   {lhs_a:10.4f}  {'<' if form_a else '>='} {m:4d}    "
              f"{lhs_b:10.6f}  {'> ' if form_b else '<='} {rhs_b:10.6f}    {s_str:>5s}")

    # Find the answer: smallest m where inequality holds for ALL k >= m
    for m in range(2, 10000):
        ln_m = math.log(m)
        if 2 * ln_m ** 2 >= m:
            first_always = None
        elif first_always is None:
            first_always = m
            # Verify it holds for the next 1000 integers
            verified = True
            for k in range(m, m + 1000):
                if 2 * math.log(k) ** 2 >= k:
                    verified = False
                    first_always = None
                    break
            if verified:
                break

    # --- Answer ---
    print("\n" + "=" * 60)
    print("ANSWER")
    print("=" * 60)
    print(f"\nBoth forms are equivalent for m >= 2:  {'YES' if all_same else 'NO'}")
    print(f"\nSmallest integer m from which 2*[ln(m)]^2 < m ALWAYS holds:")
    print(f"\n  >>> m = {first_always}")
    print(f"\n  Verification: 2*[ln({first_always})]^2 = {2 * math.log(first_always)**2:.6f} < {first_always}")
    print(f"  Previous:     2*[ln({first_always-1})]^2 = {2 * math.log(first_always-1)**2:.6f} >= {first_always-1}")
    print("=" * 60)
