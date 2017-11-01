# VBCMGH2017
Included in this repository are two files, a strict R script and a shiny web app.
The shiny web app can be found at https://vbc-standardization.shinyapps.io/vibriocidal_standardization/. 
This app can run one vibriocidal plate/.csv at a time. 
Instructions can be found on the aforementioned link.

The strict R script can run potentially unlimited vibriocidal plates/.csvs in one folder.

These programs were developed by Evan Offord and Raquel Becker, employees of Massachusetts General Hospital, Jackson 520, PI Dr. Jason B. Harris.
It implements a novel vibriocidal standardization method, intended to reduce vibriocidal variance. Publication Pending.
The source code is written in R-Project/HTML, implementing the CRAN packages rconnect, shiny, and nplr. 
All standardized values are given in logarithmic base 2.  If you have any questions or concerns, please email EOFFORD@partners.org
rlbecker@partners.org, or call 617-724-7529.

Vibriocidal Antibody Assay with a Monoclonal Antibody Standard

1.	Culture Vibrio cholerae on blood agar plates at 37OC overnight.

2.	Inoculate a loopful of bacteria from the plates in 15 ml BHI-medium in a conical flask with cotton plug.  
Incubate on a shaker at 220rpm at 37OC for 3 hours.

3.	Centrifuge the culture at 3000 rpm for 10 minutes. Decant supernatant and resuspend the pellet in 8mL sterile normal saline. 
Repeat one more time. 

4.	Adjust bacterial concentration by spectrophotometer at 600 nm to 0.3:

 		a.) Take a 1mL aliquot of your resuspended culture to make the follow dilutions:
 	 	  	i. 1:10 (100v.c + 900 saline)
 	 	 	ii. 1:15 (67vc +933 saline)
 	 	 	iii. 1:20 (50 v.c. + 950 saline)
 	 		iv. 1:25 (40 v.c. + 960 saline)
 	 		v. 1:30 (33 v.c. + 967 saline)
 		b.) Measure the OD600 in the spectrophotometer. The dilution closest to an O.D.
 	 	 	of 0.3 is your optimal dilution.
 		c.) You will need at least 150ul of optimally diluted bacteria per plate.

5.	Dilute monoclonal antibody IgG P1-A2 to desired concentration (2500 ng/mL is standard) in N. Saline 

6.	Heat-inactivate 15ul of test sera at 56OC in a water bath for 30 minutes and dilute to 1:10 (15ul sera +135ul N. Saline).

7.	Perform 2-fold serial dilutions of sample in a flat-bottom microtitre plate as follows:

a)	Dispense 25 µl of cold saline in all wells except column #2.

b)	Dispense 50 µl 1:10 diluted (15ul sera + 135ul saline) test sera in column #2. Serially dilute the sera 2-fold by 
mixing the solution in column #2, aspirating 25 µl and dispensing and mixing the sample in column #3 and so on, until column #12 
(1:10240 dilution).  Discard the last 25 µl from the last well on each row.  Keep the plates at 4oC until used.

8.	Prepare indicator (bacteria-complement-saline mixture).  The composition for each plate is as follows:
 		 Saline	           Diluted Bacteria    Complement
 		2.55 ml	      150 µl		 300 µl 	Full plate
 	 	2.125ml	      125 µl	 	 250 µl	 	6 rows (No edge affect)

Use immediately after preparation.

9.	Add 25 µl of the indicator to all wells except saline control wells. Incubate the plate on a shaker at 37OC for 1 hour (50 rev/min).

10.	Add 150 µl BHI/well.  Incubate for another 2-4 hours at 37oC without shaking until absorbance OD595 for growth control wells 
is between 0.20 to 0.28.

11.	Vibriocidal antibody titre is defined as the reciprocal of the highest serum dilutions resulting in greater 
than 50% OD reduction when compared to control wells without serum.  An increase of titre by 4-fold between acute and 
convalescent sera is considered to be significant.


Complement (Guinea pig) 				1:10 dilution within indicator

REAGENTS

Normal Saline: 9.0 g / L NaCl 

Blood Agar Plate: TSA + 5% sheep blood (Thermo Scientific, Remel; R01200)

Brain Heart Infusion Broth: Brain Heart Infusion dehydrated media (DIFCO, 237400)

Inter assay variation: Titer of batch pulled sera with batch of guinea pig complemented should not vary more
  	than two fold. 

Intra assay variation: Should be similar in duplicate rows. Maximum variation is 2-fold.


For the statistical analysis to work correctly, please follow the plate set-up listed below (on a standard 96-well plate)

Growth Control/Sample 1 (1:10)/ -> 1:2 dilution to 1:10240
Growth Control/Sample 1 (1:10)/ -> 1:2 dilution to 1:10240
Growth Control/Sample 2 (1:10)/ -> 1:2 dilution to 1:10240
Growth Control/Sample 2 (1:10)/ -> 1:2 dilution to 1:10240
Saline Control/Sample 3 (1:10)/ -> 1:2 dilution to 1:10240
Saline Control/Sample 3 (1:10)/ -> 1:2 dilution to 1:10240
Saline Control/mAb (ng/mL)/1:2 dilution to well 12G
Saline Control/mAb (ng/mL)/1:2 dilution to well 12H


