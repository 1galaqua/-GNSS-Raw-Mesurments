# Assignment 0

This repository provide a complete code to use the GNSS RAW measurements .
The repository contain 3 folder of RAW measurements :
- Driving
- Fixed
- Walking

We will use in our code the `.txt` files that contain all the headers we need.

We also have in this repository the folder gnssutils  that contain the code from johnsonmitchelld that provide a framework to use the data.

 Dependencies :
- pandas
- numpy
- navpy
- simplekml

---
To run the assignment , run in the folder -GNSS-Raw_mesurments:
```bash
python gnss_parser.py
```

After the run we obtain all the csv files required and the localisation with the coordinates from the GNSS RAW data
