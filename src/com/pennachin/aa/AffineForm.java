package com.pennachin.aa;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Represents an affine form. All symbol values are shared between all affine
 * forms.
 * 
 * @author Cassio Pennachin
 */
public class AffineForm {
	// Different form types
	public enum Type { SCALAR, FINITE, INFINITE };
	
	// Index of last noise term used
	static int lastEpsilon = 0;
	
	// Min and max kept for interval operations. r is shared noise term
	// for computation-introduced errors
	double central = 0.0;
	double min = 0.0;
	double max = 0.0;
	double r = 0.0;
	
	boolean overriden = false;
	
	Type type;
	
	HashMap<Integer, Double> terms = new HashMap<Integer, Double>();

	/**
	 * Creates a scalar form equivalent to the floating point number c.
	 */
	public AffineForm(double c) {
		if (Double.isNaN(c) || Double.isInfinite(c))
			setInfinity();
		else {
			central = c;
			min = c;
			max = c;
			type = Type.SCALAR;
		}
	}
	
	private void setInfinity() {
		type = Type.INFINITE;
		min = Double.NEGATIVE_INFINITY;
		max = Double.POSITIVE_INFINITY;
		central = 0.0;
	}

	/**
	 * Creates an affine form with the given min and max values. Checks for
	 * infinite and empty forms. Uses a single noise symbol (use -1 to 
	 * automatically allocate a new one).
	 */
	public AffineForm(double min, double max, int symbol) {
		if (Double.compare(min, max) > 0) {
			throw new RuntimeException("Trying to create affine form with inverted range: [" + min + ", " + max + "]");
		} else if (Double.isNaN(min) || Double.isNaN(max) || Double.isInfinite(min) || Double.isInfinite(max)) {
			setInfinity();
		} else {
			this.min = min;
			this.max = max;
			if (Double.compare(min, max) == 0) {
				central = min;
				type = Type.SCALAR;
			} else {
				central = (max + min) / 2.0;
				type = Type.FINITE;
				if (symbol == -1)
					symbol = ++lastEpsilon;
				terms.put(symbol, Math.max(max - central, central - min));
			}
		}
	}
	
	/**
	 * Creates an affine form from a list of terms and a given central value.
	 * r represents previous error and can be zero but shouldn't be negative.
	 */
	public AffineForm(double c, HashMap<Integer, Double> ts, double r) {
		central = c;
		this.r = r;
		
		if (Double.isNaN(c) || Double.isInfinite(c) || Double.isNaN(r) || Double.isInfinite(r)) {
			setInfinity();
			return;
		} else if (Double.compare(r, 0.0) < 0) 
			throw new RuntimeException("Trying to create affine form with negative noise: " + r);
		
		terms = new HashMap<Integer, Double>(ts);
		for (Double val : terms.values()) {
			if (Double.isNaN(val) || Double.isInfinite(val)) {
				setInfinity();
				return;
			}
		}	
		
		double radius = getRadius();
		if (Double.compare(radius, 0.0) > 0 || Double.compare(r, 0.0) > 0) {
			type = Type.FINITE;
			min = central - r - radius;
			max = central + r + radius;
		} else {
			type = Type.SCALAR;
			min = max = central;
		}
	}

	/**
	 * Creates an affine form from a list of terms and a given central value.
	 * r represents previous error and can be zero but shouldn't be negative.
	 * Enforces the provided min and max values even if they're different from
	 * the radius, allowing hybrid IA/AA model.
	 */
	public AffineForm(double c, HashMap<Integer, Double> ts, double r, double min, double max) {
		this(c, ts, r);
//		System.out.println(this);
		
		if (Double.isNaN(min) || Double.isInfinite(min) || Double.isNaN(max) || Double.isInfinite(max)) {
				setInfinity();
				return;
		}
		
		if (Double.compare(min, max) > 0) 
			throw new RuntimeException("Trying to create affine form with inverted range: [" + min + ", " + max + "]");
		
		if (Double.compare(min, this.min) > 0) {
			this.min = min;
			overriden = true;
		}
		if (Double.compare(max, this.max) < 0) {
			this.max = max;
			overriden = true;
		}
		if (Double.compare(min, max) > 0) 
			throw new RuntimeException("Trying to create affine form with inverted range: [" + min + ", " + max + "]");
		
		if (Double.compare(min, max) == 0 && type != Type.INFINITE)
			type = Type.SCALAR;
	}

	public double getCentral() {
		return central;
	}

	public double getMin() {
		return min;
	}

	public double getMax() {
		return max;
	}

	public double getR() {
		return r;
	}
	
	public Set<Integer> getNoiseSymbols() {
		return terms.keySet();
	}

	/** 
	 * Doesn't consider r and assumes terms are finite. Also doesn't consider 
	 * artificial ranges in hybrid forms.
	 */
	public double getRadius() {
		double radius = 0.0;
		for (double v : terms.values()) {
			radius += Math.abs(v);
		}
		return radius;
	}
	
	public Type getType() {
		return type;
	}

	public boolean isStrictlyPositive() {
		return Double.compare(min, 0.0) > 0; 
	}

	public boolean isStrictlyNegative() {
		return Double.compare(max, 0.0) < 0; 
	}
	
	public boolean isWeaklyPositive() {
		return Double.compare(min, 0.0) >= 0; 
	}
	
	public boolean isWeaklyNegative() {
		return Double.compare(max, 0.0) <= 0; 
	}
	
	public boolean isOverriden() {
		return overriden;
	}

	public static AffineForm copy(AffineForm af) {
		switch (af.type) {
		case SCALAR:
			return new AffineForm(af.central);
		case INFINITE:
			return new AffineForm(Double.POSITIVE_INFINITY);
		default:
			return new AffineForm(af.central, af.terms, af.r, af.min, af.max);
		}
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AffineForm other = (AffineForm) obj;
		if (type != other.type)
			return false;
		if (type == Type.SCALAR)
			return central == other.central;
		else if (type == Type.FINITE)
			return (Double.compare(central, other.central) == 0 && 
					Double.compare(r, other.r) == 0 &&
					terms.equals(other.terms)) &&
					((!overriden && !other.overriden) || 
					 (Double.compare(min, other.min) == 0 && Double.compare(max, other.max) == 0));
		else return true;
	}
	
	/** Adds two affine forms together */
	public AffineForm add(AffineForm other) {
		return linearComb(other, 1.0, 1.0, 0.0, 0.0);
	}
	
	/** Adds a (possibly negative) scalar to an affine form */
	public AffineForm scalarAdd(double delta) {
		return affine(1.0, delta, 0.0);
	}
	
	/** Multiplies an affine form by a given scalar */
	public AffineForm scalarMult(double alpha) {
		return affine(alpha, 0.0, 0.0);
	}
	
	/** Negation */
	public AffineForm negate() {
		return affine(-1.0, 0.0, 0.0);
	}
	
	/** Subtraction */
	public AffineForm sub(AffineForm other) {
		return linearComb(other, 1.0, -1.0, 0.0, 0.0);
	}
	
	/** 
	 * Multiplication. Uses the simpler approximation proposed by Stolfi et
	 * all instead of the more precise but costlier version. Computes interval
	 * product as well and keeps the intersection of the results, minimizing
	 * error propagation.
	 */
	public AffineForm mult(AffineForm other) {
		if (type == Type.INFINITE)
			return this;
		else if (other.type == Type.INFINITE)
			return other;
		else if (type == Type.SCALAR && other.type == Type.SCALAR)
			return new AffineForm(central * other.central);
		
		double c = central * other.central;
		double noise = Math.abs(central) * other.r + Math.abs(other.central) * r +
					   (getRadius() + r) * (other.getRadius() + other.r);
		
		HashMap<Integer, Double> nts = new HashMap<Integer, Double>();
		Set<Integer> symbols = new HashSet<Integer>(terms.keySet());
		symbols.addAll(other.terms.keySet());
		for (Integer sym : symbols) {
			double v1 = terms.containsKey(sym) ? terms.get(sym) : 0.0;
			double v2 = other.terms.containsKey(sym) ? other.terms.get(sym) : 0.0;
			nts.put(sym, v1 * other.central + v2 * central);
		}

		// IA multiplication for hybrid model
		List<Double> iaMult = Arrays.asList(min * other.min, min * other.max, 
											max * other.min, max * other.max);
//		System.out.println(c);
//		System.out.println(getRadius());
//		System.out.println(other.getRadius());
//		System.out.println(noise);
//		System.out.println(nts);
//		System.out.println(Collections.min(iaMult)); 
//		System.out.println(Collections.max(iaMult));
//		System.out.println("--");

		//return new AffineForm(c, nts, noise);
		return new AffineForm(c, nts, noise, Collections.min(iaMult), Collections.max(iaMult));
	}
	
	/** 
	 * Linear Combination. Returns alpha * this + beta * other + delta 
	 */
	public AffineForm linearComb(AffineForm other, double alpha, double beta, double delta, double noise) {
		if (type == Type.INFINITE)
			return this;
		else if (other.type == Type.INFINITE)
			return other;
		
		double nc = alpha * central + beta * other.central + delta;
		double nr = r + other.r + noise;
		HashMap<Integer, Double> nts = new HashMap<Integer, Double>();
		
		Set<Integer> symbols = new HashSet<Integer>(terms.keySet());
		symbols.addAll(other.terms.keySet());
		for (Integer sym : symbols) {
			double v1 = terms.containsKey(sym) ? terms.get(sym) : 0.0;
			double v2 = other.terms.containsKey(sym) ? other.terms.get(sym) : 0.0;
			nts.put(sym, alpha * v1 + beta * v2);
		}
		
		if (overriden || other.overriden) {
			double nMin = min * alpha + other.min * beta + delta;
			double nMax = max * alpha + other.max * beta + delta;
			return new AffineForm(nc, nts, nr, Math.min(nMin - noise, nMax - noise), Math.max(nMin + noise, nMax + noise));
		} else {
			return new AffineForm(nc, nts, nr);
		}
	}
	
	/**
	 * Scalar addition, multiplication and noise increment on a single form
	 */
	public AffineForm affine(double alpha, double delta, double noise) {
		if (type == Type.INFINITE)
			return this;

		double nc = central * alpha + delta;
		double nr = r * Math.abs(alpha) + noise;  
		
		HashMap<Integer, Double> nts = new HashMap<Integer, Double>(terms);
		for (Integer sym : terms.keySet()) {
			nts.put(sym, terms.get(sym) * alpha);
		}

		if (overriden) {
			double nMin = min * alpha + delta;
			double nMax = max * alpha + delta;
			return new AffineForm(nc, nts, nr, Math.min(nMin - noise, nMax - noise), Math.max(nMin + noise, nMax + noise));
		} else {
			return new AffineForm(nc, nts, nr);
		}
	}
	
	/**
	 * Exponentiation
	 */
	public AffineForm exp() {
		if (type == Type.INFINITE)
			return this;

		double iaMin = Math.exp(min);
		double iaMax = Math.exp(max);
		double delta = (iaMax + iaMin * (1.0 - min - max)) / 2.0;
		double noise = (iaMax + iaMin * (min - max - 1.0)) / 2.0;
		
		if (noise < 0 || type == Type.SCALAR) {
			double nc = Math.max(Math.exp(central), Double.MIN_VALUE);
			return new AffineForm(nc);
		}
		
		AffineForm aux = this.affine(iaMin, delta, noise);
		
		if (Double.compare(aux.min, iaMin) > 0) {
			double d = aux.min - iaMin;
			// NOTE: PLOP uses central + d, but I think that's a typo/bug, as 
			// we decrease min, so it doesn' make sense to increase central.
			return new AffineForm(aux.central - d, aux.terms, aux.r + d, iaMin, aux.max);
		} else if (Double.compare(aux.min, 0) < 0) {
			double d = Double.MIN_VALUE - aux.min;
			return new AffineForm(aux.central + d, aux.terms, aux.r + d, Double.MIN_VALUE, aux.max);
		}
		
		return aux;
	}
	
	/**
	 * Square root. Use the property that sqrt(x) = e ^ (0.5 log x)
	 * TODO: Replace with Stolfi's min-range approximation, a bit more accurate.
	 */
	public AffineForm sqrt() {
		if (type == Type.INFINITE)
			return this;
		
		AffineForm halfLog = this.log().scalarMult(0.5);
		return halfLog.exp();
	}
	
	/**
	 * Natural logarithm
	 */
	public AffineForm log() {
		if (type == Type.INFINITE)
			return this;
		// Could raise exception as well...
		else if (Double.compare(min, 0) < 0)
			return new AffineForm(Double.NEGATIVE_INFINITY);
		else if (type == Type.SCALAR)
			return new AffineForm(Math.log(central));

		double l = Math.log(min);
		double u = Math.log(max);
		double alpha = (u - l) / (max - min);
		double xs = 1 / alpha;
		double ys = (xs - min) * alpha + l;
		double logxs = Math.log(xs);
		double delta = (logxs + ys) / 2 - alpha * xs;
		double noise = Math.abs(logxs - ys) / 2;
		
		return this.affine(alpha, delta, noise);
	}
	
	/**
	 * Reciprocal, which gives us division.
	 */
	public AffineForm inv() {
		if (type == Type.INFINITE)
			return this;
		else if (type == Type.SCALAR) {
			if (Double.compare(central, 0.0) == 0)
				return new AffineForm(Double.POSITIVE_INFINITY);
			else
				return new AffineForm(1.0 / central);
		} 
		// An interval that straddles zero has ill-defined inv()
		else if (Double.compare(min, 0.0) < 0 && Double.compare(max, 0.0) > 0.0)
			return new AffineForm(Double.POSITIVE_INFINITY);
		
		double l = Math.min(Math.abs(min), Math.abs(max));
		double u = Math.max(Math.abs(min), Math.abs(max));
		double alpha = -1.0 / (u * u);
		double auxLow = 2.0 / u;
		double auxUpp = 1.0 / l - alpha * l;
		double den = Double.compare(min, 0) < 0 ? -2.0 : 2.0;
		double delta = (auxUpp + auxLow) / den;
		double noise = (auxUpp - auxLow) / 2;
		
		return this.affine(alpha, delta, Math.max(0.0, noise));		
	}
	
	/**
	 * We do division by multiplying by inv(other) as suggested by Stolfi.
	 * This means division by zero returns infinity, instead of an error.
	 */
	public AffineForm div(AffineForm other) {
		if (type == Type.INFINITE)
			return this;
		else if (other.type == Type.INFINITE)
			return other;
		
		return mult(other.inv());
	}
	
	/**
	 * Just use multiplication for the time being.
	 */
	public AffineForm sqr() {
		if (type == Type.INFINITE)
			return this;
		else if (type == Type.SCALAR)
			return new AffineForm(central * central);
		
		// Affine forms are supposed to be immutable, so we need to create a
		// new one after messing up with aux's internals.
//		System.out.println(this);
		AffineForm aux = mult(this);
//		System.out.println(aux);
		double d = aux.r - aux.central;
		if (Double.compare(d, 0.0) > 0) {
			d /= 2.0;
			aux.r -= d;
			aux.central += d;
		}
		if (Double.compare(aux.max, 0.0) > 0 && Double.compare(aux.min, 0.0) < 0) {
//			System.out.println("[" + aux.min + ", " + aux.max + "]");
			aux.max = Math.max(aux.max, -aux.min);
			aux.min = 0.0;
		}

//		System.out.println("[" + aux.min + ", " + aux.max + "]");
		return new AffineForm(aux.central, aux.terms, aux.r, aux.min, aux.max);
	}
	
	// TODO: Port proper least squares approximation
	public AffineForm sin() {
		if (type == Type.INFINITE)
			return this;
		else if (type == Type.SCALAR)
			return new AffineForm(Math.sin(central));
		else
			return new AffineForm(-1.0, 1.0, -1);
	}

	// TODO: Port proper least squares approximation
	public AffineForm cos() {
		if (type == Type.INFINITE)
			return this;
		else if (type == Type.SCALAR)
			return new AffineForm(Math.cos(central));
		else
			return new AffineForm(-1.0, 1.0, -1);
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("[");
		builder.append(central);
		builder.append(", ");
		builder.append(r);
		builder.append(", ");
		builder.append(min);
		builder.append(", ");
		builder.append(max);
		builder.append(" -- ");
		builder.append(terms);
		builder.append(", ");
		builder.append(type);
		builder.append("]");
		return builder.toString();
	}

	public String getTypeStrShort() {
		if (type == Type.INFINITE)
			return "I";
		else if (type == Type.FINITE)
			return "F";
		else return "S";
//		return String.format("[%s %.3f %.3f]", type.toString().substring(0, 1), min, max);
	}

	public boolean contains(AffineForm other) {
		if (type == Type.INFINITE)
			return true;
		else if (other.type == Type.INFINITE)
			return false;
		else
			return (min < other.min) && (max > other.max);
	}
	
	/** Simple test */
	public static void main(String[] args) {
		AffineForm x = new AffineForm(-3.4267827831536475, 0.6866020895701899, 3);
		AffineForm c1 = new AffineForm(-3.9120052);
		
		AffineForm r1 = x.negate().exp();
		AffineForm r2 = c1.sqr().sqr().exp().sqr();
		AffineForm r = r1.add(r2);

		System.out.println(r1);
		System.out.println(r2);

		// [2.682687330287812E203, 14.10199220730372, 2.682687330287812E203, 2.682687330287812E203 -- {3=-1.0350989121400223}, FINITE]
		System.out.println(r);
	}
}
