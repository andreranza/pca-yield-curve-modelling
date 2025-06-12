"""
Improved Bond Pricing Class

This class provides a cleaner, more efficient implementation for pricing
coupon bonds with better separation of concerns and lazy evaluation.
"""

from typing import List, Dict, Optional


class CouponBond:
    """
    A coupon bond with annual payments.
    
    Args:
        par_value: Face value of the bond in dollars
        maturity: Years until maturity
        coupon_rate: Annual coupon rate as decimal (e.g., 0.08 for 8%)
        price: Current market price in dollars
    """
    
    def __init__(self, par_value: float, maturity: int, coupon_rate: float, price: float):
        # Validate inputs
        if par_value <= 0:
            raise ValueError("Par value must be positive")
        if maturity <= 0:
            raise ValueError("Maturity must be positive")
        if coupon_rate < 0:
            raise ValueError("Coupon rate cannot be negative")
        if price <= 0:
            raise ValueError("Price must be positive")
            
        self.par_value = par_value
        self.maturity = maturity
        self.coupon_rate = coupon_rate
        self.price = price
        
        # Cache for expensive calculations
        self._ytm: Optional[float] = None
        self._duration: Optional[float] = None
        self._cash_flows: Optional[List[float]] = None
    
    @property
    def annual_coupon(self) -> float:
        """Annual coupon payment."""
        return self.par_value * self.coupon_rate
    
    @property
    def periods(self) -> List[int]:
        """List of payment periods (1, 2, 3, ..., maturity)."""
        return list(range(1, self.maturity + 1))
    
    def present_value_at_rate(self, discount_rate: float) -> float:
        """
        Calculate present value using given discount rate.
        
        Args:
            discount_rate: Discount rate as decimal
            
        Returns:
            Present value of the bond
        """
        if discount_rate <= -1:
            raise ValueError("Discount rate must be greater than -1")
            
        pv_coupons = sum(
            self.annual_coupon / (1 + discount_rate) ** t 
            for t in range(1, self.maturity + 1)
        )
        pv_principal = self.par_value / (1 + discount_rate) ** self.maturity
        
        return pv_coupons + pv_principal
    
    def yield_to_maturity(self, tolerance: float = 1e-6, max_iterations: int = 100) -> float:
        """
        Calculate yield to maturity using bisection method.
        
        Args:
            tolerance: Convergence tolerance
            max_iterations: Maximum number of iterations
            
        Returns:
            Yield to maturity as decimal
        """
        if self._ytm is not None:
            return self._ytm
        
        # Initial bounds
        r_low, r_high = 0.0, 1.0  # 0% to 100%
        
        # Extend upper bound if needed
        while self.present_value_at_rate(r_high) > self.price:
            r_high *= 2
            if r_high > 10:  # Sanity check
                raise ValueError("Cannot find upper bound for YTM")
        
        # Bisection search
        for _ in range(max_iterations):
            r_mid = (r_low + r_high) / 2
            pv = self.present_value_at_rate(r_mid)
            
            if abs(pv - self.price) < tolerance:
                self._ytm = r_mid
                return self._ytm
            
            if pv > self.price:
                r_low = r_mid
            else:
                r_high = r_mid
        
        raise ValueError(f"YTM did not converge within {max_iterations} iterations")
    
    def coupon_cash_flows_pv(self) -> List[float]:
        """Present value of each coupon payment."""
        if self._cash_flows is not None:
            return self._cash_flows
            
        ytm = self.yield_to_maturity()
        self._cash_flows = [
            self.annual_coupon / (1 + ytm) ** t 
            for t in range(1, self.maturity + 1)
        ]
        return self._cash_flows
    
    def duration(self) -> float:
        """Macaulay duration of the bond."""
        if self._duration is not None:
            return self._duration
            
        cash_flows = self.coupon_cash_flows_pv()
        ytm = self.yield_to_maturity()
        
        # Include present value of principal repayment
        pv_principal = self.par_value / (1 + ytm) ** self.maturity
        total_pv = sum(cash_flows) + pv_principal
        
        # Weight each cash flow by its time period
        weighted_time = sum(
            t * cf for t, cf in zip(range(1, self.maturity + 1), cash_flows)
        )
        # Add weighted time for principal
        weighted_time += self.maturity * pv_principal
        
        self._duration = weighted_time / total_pv
        return self._duration
    
    def modified_duration(self) -> float:
        """Modified duration (price sensitivity to yield changes)."""
        ytm = self.yield_to_maturity()
        return self.duration() / (1 + ytm)
    
    def total_coupon_payments(self) -> float:
        """Total coupon payments over the life of the bond."""
        return self.maturity * self.annual_coupon
    
    def interest_on_interest(self) -> float:
        """
        Interest earned by reinvesting coupons at the YTM rate.
        
        This assumes each coupon is reinvested at the YTM until maturity.
        """
        ytm = self.yield_to_maturity()
        
        # Future value of annuity formula
        if ytm == 0:
            future_value_coupons = self.total_coupon_payments()
        else:
            future_value_coupons = self.annual_coupon * (
                ((1 + ytm) ** self.maturity - 1) / ytm
            )
        
        return future_value_coupons - self.total_coupon_payments()
    
    def capital_gain(self) -> float:
        """Capital gain/loss (positive if bought below par)."""
        return self.par_value - self.price
    
    def current_yield(self) -> float:
        """Current yield (annual coupon / current price)."""
        return self.annual_coupon / self.price
    
    def summary(self) -> Dict[str, float]:
        """Return a summary of all bond metrics."""
        return {
            'par_value': self.par_value,
            'price': self.price,
            'annual_coupon': self.annual_coupon,
            'maturity': self.maturity,
            'ytm': self.yield_to_maturity(),
            'current_yield': self.current_yield(),
            'duration': self.duration(),
            'modified_duration': self.modified_duration(),
            'total_coupon_payments': self.total_coupon_payments(),
            'interest_on_interest': self.interest_on_interest(),
            'capital_gain': self.capital_gain()
        }
    
    def __repr__(self) -> str:
        return (f"CouponBond(par_value={self.par_value}, maturity={self.maturity}, "
                f"coupon_rate={self.coupon_rate:.1%}, price={self.price})")


# Example usage and testing
if __name__ == '__main__':
    # Test the examples from the original code
    print("=== Bond 1: Par Bond ===")
    bond1 = CouponBond(par_value=1000, maturity=30, coupon_rate=0.08, price=1000)
    print(f"Bond: {bond1}")
    print(f"YTM: {bond1.yield_to_maturity():.4f} ({bond1.yield_to_maturity():.2%})")
    print(f"Duration: {bond1.duration():.2f} years")
    print(f"Current Yield: {bond1.current_yield():.2%}")
    print()
    
    print("=== Bond 2: Discount Bond ===")
    bond2 = CouponBond(par_value=1000, maturity=15, coupon_rate=0.07, price=769.40)
    print(f"Bond: {bond2}")
    print(f"YTM: {bond2.yield_to_maturity():.4f} ({bond2.yield_to_maturity():.2%})")
    print(f"Duration: {bond2.duration():.2f} years")
    print(f"Capital Gain: ${bond2.capital_gain():.2f}")
    print()
    
    print("=== Bond 3: Discount Bond ===")
    bond3 = CouponBond(par_value=1000, maturity=30, coupon_rate=0.08, price=850)
    print(f"Bond: {bond3}")
    print(f"Summary: {bond3.summary()}")
    print()
    
    # Yield curve analysis
    print("=== Yield by Maturity Analysis ===")
    yields_by_maturity = []
    for mat in range(1, 31):
        bond = CouponBond(par_value=1000, maturity=mat, coupon_rate=0.08, price=950)
        yields_by_maturity.append(bond.yield_to_maturity())
    
    print("Maturity (years) -> YTM")
    for i, ytm in enumerate(yields_by_maturity[::5], 1):  # Every 5th year
        print(f"{i*5:2d} -> {ytm:.4f} ({ytm:.2%})")
