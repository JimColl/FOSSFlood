# **FOSSFlood: Free Open Source Software for mapping flood impacts anywhere**  

<table border="0">
 <tr>
    <td><img src="https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSSplash.PNG"  alt="1"></td>
    <td><img src="https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSFlood_Init.gif"  alt="1"></td>
 </tr>
</table>

**About:**  

FOSSFlood is a Free & Open Source, GUI driven application built on a combination of R, Leaflet, and Shiny that removes the barriers of entry and allows anyone to take National Water Model forecasts and transform them into actionable intelligence.  Simply...

1. download this repository
2. unzip it
3. double click the RunMe.hta, 
4. specify an area (via zip code or other common flood units), and 
5. explore the outputs!

See the [ESRI Story Maps formatted user guide](https://arcg.is/1HGCnL) for instructions on how to use the FOSSFlood application and other documentation.  

This app is fully self contained, meaning the program will only add files under the FOSSFlood-master folder. To "uninstall" FOSSFlood, just delete the folder. Pandoc and PhantomJS are also installed, and may be removed through the **add or remove programs** system dialog and the _user/AppData/Roaming_ folder respectively.

**FOSSFlood on systems other than Windows:**
**TODO**

Find or encounter an error?  Perform the following steps:
1) Close the program and delete the contents of the /AOI/ folder
1A) Ensure you have requested a valid area and entered the zip codes in the correct format
2) Restart the PC at least twice.
3) Attempt to run the program on a virtual machine
4) If the error persists, **without closing the RunMe.hta dialog**, copy the ShinyApp.log file and email it to me and we'll debug it as best we can.

**Known bugs/To Do:**  

Global file
* OpenAddress database cleanup is aggressive
* Add backups for base data download and flood mapping (archived code)

application
* shut down shiny server on idle fix
* UI flickering

printing
* Title page references featurefile (does not exist with aoi calls)

