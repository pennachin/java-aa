package com.pennachin.aa;

import static org.junit.Assert.*;

import java.util.HashMap;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

public class AffineFormTest {
	private static final double PRECISION = 0.000001;
	private AffineForm af1, af2, scl, lgr, rst;
	
	@Before
	public void setUp() {
		af1 = new AffineForm(1.0, 2.0, -1);
		af2 = new AffineForm(1.0, 2.0, 2);
		
		// Larger value
		HashMap<Integer, Double> terms = new HashMap<Integer, Double>();
		terms.put(1, 2.0);
		terms.put(2, 1.0);
		lgr = new AffineForm(10.0, terms, 0.0);
		
		// Scalar form
		scl = new AffineForm(1.0);
		
		// Restricted range through manual min/max values
		terms = new HashMap<Integer, Double>();
		terms.put(2, 0.5);
		rst = new AffineForm(1.5, terms, 0.0, 1.1, 1.9);
	}
	
	@Test
	public void testAddForms() {
		AffineForm sum = af1.add(af2);
		
		Set<Integer> noiseSymbols = sum.getNoiseSymbols();
		double termValueSum = 0.0;
		for (Double v : sum.terms.values()) {
			termValueSum += v;
		}
		
		assertEquals(3.0, sum.getCentral(), PRECISION);
		assertEquals(2.0, sum.getMin(), PRECISION);
		assertEquals(4.0, sum.getMax(), PRECISION);
		assertEquals(0.0, sum.getR(), PRECISION);
		assertEquals(1.0, sum.getRadius(), PRECISION);
		
		assertEquals(2, noiseSymbols.size());
		assertTrue(noiseSymbols.contains(1));
		assertTrue(noiseSymbols.contains(2));
		assertEquals(termValueSum, sum.getRadius(), PRECISION);
	}
	
	@Test
	public void testAddFormsCustomRange() {
		AffineForm sum = af1.add(rst);
		
		assertEquals(3.0, sum.getCentral(), PRECISION);
		assertEquals(2.1, sum.getMin(), PRECISION);
		assertEquals(3.9, sum.getMax(), PRECISION);
		assertEquals(0.0, sum.getR(), PRECISION);
		// Doesn't change with artificial range!
		assertEquals(1.0, sum.getRadius(), PRECISION);
	}

	@Test
	public void testAddScalars() {
		AffineForm sum1 = af1.scalarAdd(1.0);
		AffineForm sum2 = af1.add(scl);
		
		assertEquals(2.5, sum1.getCentral(), PRECISION);
		assertEquals(2.0, sum1.getMin(), PRECISION);
		assertEquals(3.0, sum1.getMax(), PRECISION);
		assertEquals(0.0, sum1.getR(), PRECISION);
		assertEquals(0.5, sum1.getRadius(), PRECISION);
		
		assertTrue(sum1.equals(sum2));
	}
	
	@Test
	public void testAddInfinite() {
		AffineForm af3 = new AffineForm(Double.POSITIVE_INFINITY);
		AffineForm sum1 = af3.add(af1);
		AffineForm sum2 = af1.add(af3);
		
		assertEquals(AffineForm.Type.INFINITE, af3.type);
		assertEquals(AffineForm.Type.INFINITE, sum1.type);
		assertEquals(AffineForm.Type.INFINITE, sum2.type);
		assertEquals(sum1, sum2);
	}

	@Test
	public void testMultScalar() {
		AffineForm prod1 = af1.scalarMult(0.5);
		AffineForm prod2 = af1.scalarMult(Double.NaN);
		AffineForm prod3 = scl.scalarMult(0.5);
		AffineForm prod4 = scl.scalarMult(Double.NaN);
		
		assertEquals(AffineForm.Type.FINITE, prod1.type);
		assertEquals(AffineForm.Type.INFINITE, prod2.type);
		assertEquals(AffineForm.Type.SCALAR, prod3.type);
		assertEquals(AffineForm.Type.INFINITE, prod4.type);
		
		assertEquals(0.75, prod1.getCentral(), PRECISION);
		assertEquals(0.25, prod1.getRadius(), PRECISION);
		assertEquals(0.5, prod3.getCentral(), PRECISION);
	}

	@Test
	public void testNegation() {
		AffineForm neg = af1.negate();
		
		assertEquals(-1.5, neg.getCentral(), PRECISION);
		assertEquals(-2.0, neg.getMin(), PRECISION);
		assertEquals(-1.0, neg.getMax(), PRECISION);
		assertEquals(0.0, neg.getR(), PRECISION);
		assertEquals(0.5, neg.getRadius(), PRECISION);		
	}
	
	@Test
	public void testLinearCombination() {
		AffineForm lc = af1.linearComb(af2, 2.0, 0.5, -10.0, 0.0);

		assertEquals(-6.25, lc.getCentral(), PRECISION);
		assertEquals(-7.5, lc.getMin(), PRECISION);
		assertEquals(-5.0, lc.getMax(), PRECISION);
		assertEquals(0.0, lc.getR(), PRECISION);
		assertEquals(1.25, lc.getRadius(), PRECISION);		
	}
	
	@Test
	public void testMultiplication() {
		HashMap<Integer, Double> terms1 = new HashMap<Integer, Double>();
		terms1.put(1, 2.0);
		terms1.put(2, 1.0);
		
		HashMap<Integer, Double> terms2 = new HashMap<Integer, Double>();
		terms2.put(1, -2.0);
		terms2.put(3, 1.0);

		AffineForm af3 = new AffineForm(10.0, terms2, 0.0);
		AffineForm mult = lgr.mult(af3);
		
		assertEquals(AffineForm.Type.FINITE, mult.type);
		assertEquals(100.0, mult.getCentral(), PRECISION);
		assertEquals(71.0, mult.getMin(), PRECISION);
		assertEquals(129.0, mult.getMax(), PRECISION);
		assertEquals(9.0, mult.getR(), PRECISION);
	}
	
	@Test
	public void testExp() {
		// Test scalar
		AffineForm af3 = new AffineForm(3.5);
		AffineForm af4 = new AffineForm(-1.0); 
		AffineForm exp1 = scl.exp();
		AffineForm exp2 = af3.exp();
		AffineForm exp3 = af4.exp();
		
		assertEquals(Math.E, exp1.getCentral(), PRECISION);
		assertEquals(Math.exp(3.5), exp2.getCentral(), PRECISION);
		assertEquals(Math.exp(-1.0), exp3.getCentral(), PRECISION);
		
		// Test interval
		AffineForm exp4 = af1.exp();
		
		assertEquals(5.06, exp4.getCentral(), 0.01);
		assertEquals(0.98, exp4.getR(), 0.01);
		assertEquals(2.72, exp4.getMin(), 0.01);
		assertEquals(7.39, exp4.getMax(), 0.01);

		AffineForm exp5 = lgr.exp();

		assertEquals(221755.0, exp5.getCentral(), 1.0);
		assertEquals(217368.0, exp5.getR(), 1.0);
		assertEquals(1097.0, exp5.getMin(), 1.0);
		assertEquals(442413.0, exp5.getMax(), 1.0);
	}
	
	@Test
	public void testLog() {
		AffineForm log1 = af1.log();
		assertEquals(0.38, log1.getCentral(), 0.01);
		assertEquals(0.02, log1.getR(), 0.01);
		assertEquals(0.0, log1.getMin(), 0.01);
		assertEquals(0.75, log1.getMax(), 0.01);
		
		AffineForm log2 = lgr.log();
		assertEquals(2.28, log2.getCentral(), 0.01);
		assertEquals(0.02, log2.getR(), 0.01);
		assertEquals(1.95, log2.getMin(), 0.01);
		assertEquals(2.61, log2.getMax(), 0.01);

		AffineForm log3 = rst.log();
		assertEquals(0.39, log3.getCentral(), 0.01);
		assertEquals(0.02, log3.getR(), 0.01);
		assertEquals(0.09, log3.getMin(), 0.01);
		assertEquals(0.68, log3.getMax(), 0.01);
	}
	
	@Test
	public void testSqrt() {
		AffineForm sqrt1 = af1.sqrt();
		assertEquals(1.23, sqrt1.getCentral(), 0.01);
		assertEquals(0.05, sqrt1.getR(), 0.01);
		assertEquals(1.0, sqrt1.getMin(), 0.01);
		assertEquals(1.45, sqrt1.getMax(), 0.01);
		
		AffineForm sqrt2 = lgr.sqrt();
		assertEquals(3.17, sqrt2.getCentral(), 0.01);
		assertEquals(0.11, sqrt2.getR(), 0.01);
		assertEquals(2.65, sqrt2.getMin(), 0.01);
		assertEquals(3.69, sqrt2.getMax(), 0.01);

		AffineForm sqrt3 = rst.sqrt();
		assertEquals(1.23, sqrt3.getCentral(), 0.01);
		assertEquals(0.04, sqrt3.getR(), 0.01);
		assertEquals(1.05, sqrt3.getMin(), 0.01);
		assertEquals(1.40, sqrt3.getMax(), 0.01);
	}
	
	@Test
	public void testInv() {
		// Around zero
		AffineForm af3 = new AffineForm(-2.0, 2.0, -1);
		AffineForm inv = af3.inv();
		assertEquals(AffineForm.Type.INFINITE, inv.getType());
		
		// Infinity should be preserved
		AffineForm inf = new AffineForm(Double.POSITIVE_INFINITY);
		AffineForm inv2 = inf.inv();
		assertEquals(AffineForm.Type.INFINITE, inv2.getType());
		
		// Regular
		AffineForm inv3 = af1.inv();
		assertEquals(0.75, inv3.getCentral(), PRECISION);
		assertEquals(0.5, inv3.getMin(), PRECISION);
		assertEquals(1.0, inv3.getMax(), PRECISION);
		assertEquals(0.125, inv3.getR(), PRECISION);
		assertEquals(0.125, inv3.getRadius(), PRECISION);
	}
	
	@Test 
	public void testDiv() {
		// Div by zero should return infinity
		AffineForm zero = new AffineForm(0.0);
		AffineForm div = af1.div(zero);
		assertEquals(AffineForm.Type.INFINITE, div.getType());
		
		// Regular division is tested by inversion + multiplication
	}
	
	@Test 
	public void testSqr() {
		AffineForm mult1 = af1.mult(af1);
		AffineForm sqr1  = af1.sqr();
		AffineForm mult2 = lgr.mult(lgr);
		AffineForm sqr2  = lgr.sqr();
		AffineForm mult3 = rst.mult(rst);
		AffineForm sqr3  = rst.sqr();
		
		assertEquals(mult1, sqr1); // small interval
		assertEquals(mult2, sqr2); // larger, multiple symbols
		assertEquals(mult3, sqr3); // artificially restricted bounds
	}
	
	// Compute g(x) = sqrt(x^2 - x + 0.5) / sqrt(x^2 + 0.5) over 16 intervals in [-2, 2]
	// TODO: Enable this with new golden standard after replacing sqrt() function.
	@Test
	public void testChainedExpressions() {
		// Center, noise, min, max for each interval.
//		double[][] truthTable = {
//				{ 1.2151, 0.0345, 1.1802, 1.2499 },
//				{ 1.2367, 0.0420, 1.1930, 1.2803 },
//				{ 1.2610, 0.0522, 1.2049, 1.3172 },
//				{ 1.2868, 0.0411, 1.2372, 1.3364 },
//				{ 1.3088, 0.0498, 1.2420, 1.3757 },
//				{ 1.3121, 0.0554, 1.2229, 1.4013 },
//				{ 1.2628, 0.0507, 1.1509, 1.3747 },
//				{ 1.1151, 0.0427, 0.9991, 1.2311 },
//				{ 0.8726, 0.0503, 0.7447, 1.0139 },
//				{ 0.6470, 0.0396, 0.5744, 0.7477 },
//				{ 0.5504, 0.0365, 0.5075, 0.5933 },
//				{ 0.5608, 0.0416, 0.4971, 0.6246 },
//				{ 0.6072, 0.0387, 0.5410, 0.6734 },
//				{ 0.6558, 0.0325, 0.5978, 0.7139 },
//				{ 0.6980, 0.0266, 0.6494, 0.7466 },
//				{ 0.7329, 0.0218, 0.6925, 0.7732 }
//		};
//		
//		for (int i = 0; i < 16; i++) {
//			double lower = -2.0 + i * 0.25;
//			double upper = lower + 0.25;
//			AffineForm input = new AffineForm(lower, upper, -1);
//			
//			AffineForm square = input.sqr();
//			AffineForm num2 = square.sub(input).scalarAdd(0.5);
//			AffineForm den2 = square.scalarAdd(0.5);
//			AffineForm num = num2.sqrt();
//			AffineForm den = den2.sqrt();
//			AffineForm result = num.div(den);
			
//			System.out.println(result);
			
//			assertEquals(truthTable[i][0], result.getCentral(), 0.0001);
//			assertEquals(truthTable[i][1], result.getR(), 0.0001);
//			assertEquals(truthTable[i][2], result.getMin(), 0.0001);
//			assertEquals(truthTable[i][3], result.getMax(), 0.0001);
//		}
	}
	
	@Test
	public void testCreation() {
		assertEquals(af1.getCentral(), af2.getCentral(), PRECISION);
		assertEquals(af1.getMin(), af2.getMin(), PRECISION);
		assertEquals(af1.getMax(), af2.getMax(), PRECISION);
		assertEquals(af1.getR(), af2.getR(), PRECISION);
		assertEquals(af1.getRadius(), af2.getRadius(), PRECISION);

		assertEquals(1.5, af2.getCentral(), PRECISION);
		assertEquals(1.0, af2.getMin(), PRECISION);
		assertEquals(2.0, af2.getMax(), PRECISION);
		assertEquals(0.0, af2.getR(), PRECISION);
		assertEquals(0.5, af2.getRadius(), PRECISION);
		
		assertEquals(1.0, scl.getCentral(), PRECISION);
		assertEquals(scl.getMin(), scl.getMax(), PRECISION);
		assertEquals(0.0, scl.getRadius(), PRECISION);
	}
	
	@Test
	public void testEquality() {
		AffineForm af3 = new AffineForm(1.0, 2.0, 2);
		AffineForm af4 = AffineForm.copy(af2);
		
		assertFalse(af1.equals(af2));
		assertTrue(af3.equals(af2));		
		assertTrue(af4.equals(af2));		
	}
	
	@Test 
	public void testAroundZero() {
		AffineForm af3 = new AffineForm(-1.0, 2.0, -1);
		AffineForm af4 = new AffineForm(0.0, 2.0, -1);
		AffineForm af5 = new AffineForm(0.0000001, 2.0, -1);
		AffineForm af6 = new AffineForm(-2.0, 0.0, -1);
		AffineForm af7 = new AffineForm(-2.0, -0.0000001, -1);
		
		assertFalse(af3.isWeaklyPositive());
		assertFalse(af3.isWeaklyNegative());
		
		assertTrue(af4.isWeaklyPositive());
		assertFalse(af4.isStrictlyPositive());
		assertFalse(af4.isWeaklyNegative());
		
		assertTrue(af5.isWeaklyPositive());
		assertTrue(af5.isStrictlyPositive());
		assertFalse(af5.isWeaklyNegative());
		
		assertTrue(af6.isWeaklyNegative());
		assertFalse(af6.isStrictlyNegative());
		assertFalse(af6.isWeaklyPositive());
		
		assertTrue(af7.isWeaklyNegative());
		assertTrue(af7.isStrictlyNegative());
		assertFalse(af7.isWeaklyPositive());
	}
	
	@Test(expected = RuntimeException.class)
	public void testConstructionExceptionRange() {
		new AffineForm(2.0, 1.0, -1);
	}
	
	@Test
	public void testInfiniteTerm() {
		HashMap<Integer, Double> terms = new HashMap<Integer, Double>();
		terms.put(1, Double.POSITIVE_INFINITY);
		AffineForm af3 = new AffineForm(0.0, terms, 0.0);
		assertEquals(AffineForm.Type.INFINITE, af3.type);
	}
}
