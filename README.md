# SatelliteGeodesyInCPP

## General info
- This repository contains some of the most useful algorithms and OOP classes in Satellite Positioning.

- The projects in the repository were created while I was teaching coding in C++ to myself. 

  (To be honest, it helped me a lot to not only in practicing C++ but also reviewing my knowledge in Satellite Geodesy)

- The main idea behind the project was to create an application which performs single point positioning in C++.

- Each subfolder in the project is created for only a step in GPS/GNSS data processing.

### Details:
#### Adjustment Module
Performs Least Squares Adjustment of the GPS/GNSS data epochwise by creating the design, clsoure matrices and performs the matrix oparations with the help of Matrix class in the same repo.

#### BrdcEph
Calculates the position of any GPS/GNSS satellite at any epoch given the rinex navigation file as input.

#### Constants
Contains the constant values in the processing of GPS/GNSS data.

#### Ellipsoid
Defines an ellipsoid which is to be used in the whole project given the ellipsoid parameters.

#### GPSModule
Contains basic utilities for coordinate transformations or some other useful helping calculations such as interplation of sv clocks, troposheric corrections for processing GPS/GNSS data. And some more ...

#### RinexN
Contains utility for reading rinex navigation file and some more including estimation and creating plots...

#### RinexObsFile
Contains utility for reading rinex observation file and some more including creating plots.

#### SinglePointPositioning
Driver program for the project

#### TimeSystems
Contains time conversions in GPS/GNSS data processing.
