# **FOSSFlood: Free Open Source Software for mapping flood impacts anywhere**  

See the [**ESRI Story Maps formatted user guide**](https://arcg.is/1HGCnL) for instructions on how to use the FOSSFlood application and other documentation.  

<table border="0">
 <tr>
    <td><img src="https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSSplash.PNG"  alt="1"></td>
    <td><img src="https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSFlood_Init.gif"  alt="1"></td>
 </tr>
</table>


**About:**  
FOSSFlood is a Free & Open Source, GUI driven application built on a combination of R, Leaflet, and Shiny that removes the barriers of entry and allows anyone to take National Water Model forecasts and transform them into actionable intelligence.  

**I'm in a rush, TL;DR?**
FOSSFlood is Free and Open Source, GUI based, accessible, flood impact mapping software.  To execute it:
1. download this repository
2. unzip it
3. double click the RunMe.hta, 
4. specify an area (via zip code or other common flood units), and 
5. explore the outputs!

**Troubleshooting:**  
Find or encounter an error?  Perform the following steps:
1) Close the program and delete the contents of the /AOI/ folder.
1A) Ensure you have requested a valid area and entered the zip codes in the correct format.
2) Restart the PC at least twice.
3) Attempt to run the program on a virtual machine.
4) If the error persists, **without closing the RunMe.hta dialog**, copy the ShinyApp.log file and email it to me and we'll debug it as best we can.

**Known bugs/To Do:**  
Global file:
* NWIS cleaning, validation, error checking
* OpenAddress database cleanup is aggressive
* Add backups for base data download and flood mapping (archived code)

Application:
* shut down shiny server on idle fix
* UI flickering

Printing:
* Title page references featurefile (does not exist with aoi calls)