"""
This script semiautomate the collection of bond data from
'https://markets.businessinsider.com/bonds/finder'.

It returns a csv file having as columns the bond features.
"""

import pandas as pd

def bondsDataFetcher():

	bonds = {
	'name': list(),
	'mrkt_price': list(),
	'cp_rate': list(),
	'ytm': list(),
	'isn': list(),
	'issue_price': list(),
	'issue_date': list(),
	'face_val': list(),
	'maturity': list(),
	'coupon_pymt_date': list(),
	'numb_payments': list(),
	'startcoupon_date': list(),
	'finalcoupon_date': list()
	}


	active = True

	while active:

		name = str(input("\nBOND NAME: "))
		bonds['name'].append(name)


		try:
			price = float(input("\nMARKET PRICE ($): "))
		except ValueError:
			print("ValueError: invalid price. Please, try again.")
			price = float(input("\nMARKET PRICE: "))
		else:
			bonds['mrkt_price'].append(price)


		try:
			coupon_rate = float(input("\nCOUPON RATE (1.75%): "))
		except ValueError:
			print("ValueError: invalid coupon rate. Please, try again.")
			coupon_rate = float(input("\nCOUPON RATE (1.75%): "))
		else:
			bonds['cp_rate'].append(coupon_rate)


		try:
			ytm = float(input("\nYTM (1.57%): "))
		except ValueError:
			print("ValueError: invalid yield-to-maturity. Please, try again.")
			ytm = float(input("\nYTM (1.57%): "))
		else:
			bonds['ytm'].append(ytm)


		isn = str(input("\nISN: "))
		bonds['isn'].append(isn)


		try:
			issue_price = float(input("\nISSUE PRICE ($): "))
		except ValueError:
			print("ValueError: invalid issue price. Please, try again.")
			issue_price = float(input("\nISSUE PRICE: "))
		else:
			bonds['issue_price'].append(issue_price)


		issue_date = str(input("\nISSUE DATE (MM-DD-YYYY): "))
		bonds['issue_date'].append(issue_date)


		try:
			face = float(input("\nDENOMINATION FACEVALUE ($): "))
		except ValueError:
			print("ValueError: invalid facevalue. Please, try again.")
			face = float(input("\nDENOMINATION FACEVALUE: "))
		else:
			bonds['face_val'].append(face)


		maturity = str(input("\nMATURITY date (MM-DD-YYYY): "))
		bonds['maturity'].append(maturity)


		coupon_pymt_date = str(input("\nCOUPON PAYMENT DATE (MM-DD-YYYY): "))
		bonds['coupon_pymt_date'].append(coupon_pymt_date)


		try:
			numb_payments = int(input("\nNUMB OF PAYMs X YR: "))
		except ValueError:
			print("ValueError: ivalid number of payments. Please, try again.")
			numb_payments = int(input("\nNUMB OF PAYMs X YR: "))
		else:
			bonds['numb_payments'].append(numb_payments)

		startcoupon_date = str(input("\nCOUPON START DATE (MM-DD-YYYY): "))
		bonds['startcoupon_date'].append(startcoupon_date)

		finalcoupon_date = str(input("\nFINAL COUPON DATE (MM-DD-YYYY): "))
		bonds['finalcoupon_date'].append(finalcoupon_date) 


		message = input("\nDO U WANT TO CONTINUE? (yes/no): ")

		if message == 'no':
			active = False

	bonds_df = pd.DataFrame.from_dict(bonds)
	bonds_df.to_csv("../data/bonds-data-prova.csv")

	return bonds