Here you can find the raw data we used for our analysis.

As of 21/Dec/2023 raw data are stored in the EDI [repository](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-luq.233.2). We modified the code in order users to directly pull the data to their local machine.

## Metadata
- ```Tag``` = Unique tag number for the individual stem. Tag number is the relational variable which we can match with historic data of LDFP.

- ```basal_d``` = diameter at base of stem (cm)

- ```dbh``` = diameter at breast height (cm) 

- ```height``` = stem height (m)

- ```exposure``` = visual inspection of canopy exposure (factors : 1-4) 

- ```h_dbh``` = height at which dbh was measured (m)

- ```notes```= notes taken in the field

- ```second coord``` = second cordinate; some palms either had missing tags or were located outside the boundaries of LFDP. In case an individual was close to a known (taged) palm, then the unique tag number of the closest in proximity palm is stored here.

For each of these palms, we measured stem height from the ground to the base of the crown (H<sub>bc</sub>; height of the youngest internode), diameter at 130 cm above ground (D<sub>130</sub>), and basal diameter (D<sub>B</sub>; just above the top of the roots) (Figure S1). For palms with H<sub>bc</sub> ≤ 2 m, we measured height with a measuring tape, for palms between 2 and 5 m tall, we used a calibrated extension pole, and for palms > 5 m, we used a Nikon forestry Pro II rangefinder. ```Second coordinate``` was used to associate palms of unknown tag with environmental heterogeneity dataset. ```Exposure``` is not part of our analysis.

Update:
EDI introduced an additional column with date. We notify the users that this is just an approximate date of the sampling date.
