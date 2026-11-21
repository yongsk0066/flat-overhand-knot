"""Exact separating Farkas certificate for both pentagon taco branches.

A vector d with grad(core).d < 0 and grad(w_i).d > 0 for all eight
walls rules out grad(core)=sum(lambda_i grad(w_i)), lambda_i>=0,
on either seven-wall branch. This direct alternative avoids a numerical
search over multiplier signs. Band orientations and all signs are checked
in exact algebra / outward-rounded intervals by exact_core_descent.
"""
from exact_core_descent import certify

if __name__ == '__main__':
    certify()
