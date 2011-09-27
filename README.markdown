# Java Affine Arithmetic Library

This is a simple implementation of affine arithmetic in Java. It is inspired by two previous implementations:

* [LibAffa](http://www.nongnu.org/libaffa/), in C++
* [PLOP](http://code.google.com/p/plop), a Common Lisp probabilistic program evolution system which includes affine arithmetic for detection and rejection of infinite intervals.

Affine arithmetic is an interval method developed in 1993. The best (lengthy) description is this [volume](http://www.dcc.unicamp.br/~stolfi/EXPORT/papers/by-tag/fig-sto-97-iaaa.ps.gz).

# Example Usage

From the main function in AffineForm:

    // A finite interval
    AffineForm x = new AffineForm(-3.4267827831536475, 0.6866020895701899, 3);

    // A scalar
    AffineForm c1 = new AffineForm(-3.9120052);

    // exp(-x)
    AffineForm r1 = x.negate().exp();
    
    // Chained operations
    AffineForm r2 = c1.sqr().sqr().exp().sqr();
    
    // Operations with two arguments
    AffineForm r = r1.add(r2);

# TODO

* Proper Ant and JUnit integration. Not a Java developer by trade, I'm just using the Eclipse-generated file, with a minor correction. You need to have lib/junit.jar to run tests.
* Proper sqrt() implementation, which should be a bit more accurate than current code.
* Proper least squares approximation for sin() and cos(), based on LibAffa.