package com.pennachin.aa;

/**
 * Represents an interval for interval arithmetic operations. Includes to/from
 * AffineForm conversion.
 * 
 * @author Cassio Pennachin
 */

// TODO: Consider rounding issues in range, radius and so forth

public class Interval implements Cloneable {
	private double min, max;

	public Interval(double max, double min) {
		if (min > max) 
			throw new RuntimeException("Trying to create interval with inverted range");
		
		this.max = max;
		this.min = min;
	}
	
	public double getMin() {
		return min;
	}

	public void setMin(double min) {
		if (min > max)
			throw new RuntimeException("Trying to invert interval range");

		this.min = min;
	}

	public double getMax() {
		return max;
	}

	public void setMax(double max) {
		if (min > max)
			throw new RuntimeException("Trying to invert interval range");

		this.max = max;
	}

	// TODO: Will this work with infinities?
	public double range() {
		return max - min;
	}
	
	public boolean straddlesZero() {
		return (min <= 0.0 && max >= 0.0);
	}
	
	public boolean isFinite() {
		return !Double.isInfinite(min) && !Double.isNaN(min) && !Double.isInfinite(max) && !Double.isNaN(max);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Interval other = (Interval) obj;
		return ((Double.compare(min, other.min) == 0) && (Double.compare(max, other.max) == 0));
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
	
	public AffineForm toAffineForm() {
		return new AffineForm(min, max, -1);
	}
}