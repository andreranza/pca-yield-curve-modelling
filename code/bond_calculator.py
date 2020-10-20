"""
This a Python Class Object which allows to instanciate coupon bonds 
given the following parameters:

- face value in dollars, 
- maturity in years, 
- annual coupon rate in decimal points,
- trading price in dollars.


It provides methods to calculate the following financial quntities:

- yield to maturity (using the bisection search algorithm),
- present value at time t=0 of each coupon paid yearly,
- interest-on-interest under the assumption of yearly yield to maturity as reinvestment rate,
- duration of the bond,
- capital gain.
"""

class multipleCouponBond:
	"""
	Create a Bond with specifics expressed on annual basis.
	price: selling price on the primary market.
	par_value: nominal face value of the bond. If price = par_value, bond sold at par.  
	mat: number of years until the bond matures.
	coupon_rate: percentage in decimals of the par_value.
	price: valiue paid on the primary market to purchase the bond
	"""

	def __init__(self, par_value, mat, coupon_rate, price):
		self.par_value = par_value
		self.mat = mat 
		self.coupon_rate = coupon_rate
		self.price = price
		self.getCoupon() 
		self.ytm()
		self.cashFlows()
		self.getActualPrice()
		self.getPeriods()
		self.totalCouponPayments()
		self.duration()
		self.interestOnInterest()
		self.capitalGain()

	def getPeriods(self):
		"""Return a list with the time periods."""
		self.periods = [t for t in range(1,self.mat+1)]
		return self.periods

	def getCoupon(self):
		"""Compute the annual coupon."""
		self.coupon = self.par_value*self.coupon_rate
		return self.coupon

	def ytm(self, r_low = 0.00, r_high = 100.00, epsilon = 0.01):
		"""Compute the yield to maturity by means of bisection search."""
		r = (r_high + r_low)/2
		# Calculate the present value of coupons
		cash_flows = [(self.coupon/((1+r)**t)) for t in range(1, self.mat+1)]
		pv_c = sum(cash_flows)
		# Calculate the present value of par value
		pv_par = self.par_value/((1+r)**(self.mat))
		# Bond value 
		pv = pv_c + pv_par

		while abs(pv - self.price) > epsilon:
			if pv < self.price: # r should decrease
				r_high = r
			else:	# r should decrease
				r_low = r
			# update r
			r = (r_high + r_low)/2
			# update pv
			cash_flows = [(self.coupon/((1+r)**t)) for t in range(1, self.mat+1)]
			pv_c = sum(cash_flows)
			pv_par = self.par_value/((1+r)**(self.mat))
			pv = pv_c + pv_par

		self.results = {'ytm':r, 'pv':pv, 'pv_c':pv_c, 
						'pv_par':pv_par, 'cash_flows':cash_flows}
		
		return self.results['ytm']

	def getActualPrice(self):
		"""Return the present value the bond."""
		return round(self.results['pv'], 2)

	def actualParValue(self):
		"""Return the present value of the final payback on the capital invested."""
		return round(self.results['pv_par'], 2)

	def cashFlows(self):
		"""Return a list with the present value at time t=0 of the coupons payments."""
		cf = map(lambda x: round(x, 2), self.results['cash_flows'])
		return list(cf)

	def duration(self):
		"""Compute the duration of the bond."""
		self.d_num = [p*cf for p,cf in zip(self.periods, self.results['cash_flows'])]
		self.d = sum(self.d_num)/self.results['pv']
		return round(self.d, 2)

	def totalCouponPayments(self):
		"""Compute the total coupon payments of the bond."""
		self.tot_coup_paym = self.mat*self.coupon
		return round(self.tot_coup_paym, 2)

	def interestOnInterest(self):
		"""Compute the interest on interest of the bond."""
		self.int_on_int = ((self.coupon)*((((1+self.results['ytm'])**(self.mat)-1))/(self.results['ytm'])))-(self.tot_coup_paym)
		return round(self.int_on_int, 2)

	def capitalGain(self):
		"""Return the capital gain of the bond."""
		self.capital_gain = self.par_value - self.price
		return self.capital_gain

	
if __name__ == '__main__':

	bond1 = multipleCouponBond(par_value = 1000, mat = 30, coupon_rate = 0.08, price = 1000)
	print('ytm:', bond1.ytm())
	print('coupon value:', bond1.getCoupon())
	print('periods:', bond1.getPeriods())
	print('duration:', bond1.duration())
	print('coupon cash flows:', bond1.cashFlows())
	print('present value:', bond1.getActualPrice())
	print('present value of par value:', bond1.actualParValue())
	print('total coupon payments:', bond1.totalCouponPayments())
	print('interest on interest:', bond1.interestOnInterest())

	bond2 = multipleCouponBond(par_value = 1000, mat = 15, coupon_rate = 0.07, price = 769.40)
	print('ytm:', bond2.ytm())
	print('coupon value:', bond2.getCoupon())
	print('periods:', bond2.getPeriods())
	print('duration:', bond2.duration())
	print('coupon cash flows:', bond2.cashFlows())
	print('present value:', bond2.getActualPrice())
	print('present value of par value:', bond2.actualParValue())
	print('total coupon payments:', bond2.totalCouponPayments())
	print('interest on interest:', bond2.interestOnInterest())
	print('capital gain:', bond2.capitalGain())
	
	bond3 = multipleCouponBond(par_value = 1000, mat = 30, coupon_rate = 0.08, price = 850)
	print('ytm:', bond3.ytm())
	print('coupon value:', bond3.getCoupon())
	print('periods:', bond3.getPeriods())
	print('duration:', bond3.duration())
	print('coupon cash flows:', bond3.cashFlows())
	print('present value:', bond3.getActualPrice())
	print('present value of par value:', bond3.actualParValue())
	print('total coupon payments:', bond3.totalCouponPayments())
	print('interest on interest:', bond3.interestOnInterest())
	print('capital gain:', bond3.capitalGain())

	# we calculate the yield to maturity for a coupon bond with fixed parameters
	# but different maturities
	yields_given_mat = list()
	for mat in range(1,31):
		yields_given_mat.append(multipleCouponBond(1000,mat,0.08,950).ytm())
	print(yields_given_mat)
