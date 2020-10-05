
bonds = {
	'name' : list(),
	'mrkt_price' : list(),
	'cp_rate' : list(),
	'ytm' : list(),
	'isn' : list(),
	'issue_price' : list(),
	'issue_date': list(),
	'face_val' : list(),
	'maturity' : list(),
	'coupon_pymt_date' : list(),
	'numb_payments' : list(),
	'startcoupon_date' : list(),
	'finalcoupon_date' : list()
}


active = True

while active:
	message = input("\nDO U WANT TO CONTINUE? (yes/no): ")
	if message == 'no':
		active = False
	else:
		name = str(input("\nBOND NAME: "))
		bonds['name'].append(name)

		price = float(input("\nMARKET PRICE: "))
		bonds['mrkt_price'].append(price)

		coupon_rate = float(input("\nCOUPON RATE (1.75%): "))
		bonds['cp_rate'].append(coupon_rate)

		ytm = float(input("\nYTM (1.57%): "))
		bonds['ytm'].append(ytm)

		isn = str(input("\nISN: "))
		bonds['isn'].append(isn)

		issue_price = float(input("\nISSUE PRICE: "))
		bonds['issue_price'].append(issue_price)

		issue_date = str(input("\nISSUE DATE: "))
		bonds['issue_date'].append(issue_date)

		face = float(input("\nDENOMINATION FACEVALUE: "))
		bonds['face_val'].append(face)

		maturity = str(input("\nMATURITY date (MM-DD-YR): "))
		bonds['maturity'].append(maturity)

		coupon_pymt_date = str(input("\nCOUPON PAYMENT DATE: "))
		bonds['coupon_pymt_date'].append(coupon_pymt_date)

		numb_payments = str(input("\nNUMB OF PAYMs X YR: "))
		bonds['numb_payments'].append(numb_payments)

		startcoupon_date = str(input("\nCOUPON START DATE (MM-DD-YR): "))
		bonds['startcoupon_date'].append(startcoupon_date)

		finalcoupon_date = str(input("\nFINAL COUPON DATE (MM-DD-YR): "))
		bonds['finalcoupon_date'].append(finalcoupon_date) 


print(bonds)

file = open("bonds_data","w")
file.write(str(bonds))
file.close()

import pandas as pd
with open("bonds_data.txt") as bonds:
    data = bonds.read()    
data = eval(data)
data = pd.DataFrame(data)

bonds.df.rename(columns=str.upper,inplace = True)
bonds.df
