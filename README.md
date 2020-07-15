# **FOSSFlood: Free Open Source Software for mapping flood impacts anywhere**

![fossflood banner](https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSSplash.PNG)
![fossflood gif](https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSFlood_Init.gif)

<table border="0">
 <tr>
    <td>![fossflood banner](https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSSplash.PNG)</td>
    <td>![fossflood gif](https://github.com/JimColl/FOSSFlood/blob/master/data/misc/FOSSFlood_Init.gif)</td>
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

**Known bugs/To Do:**  

Global file
* AOI's crossing state lines are broken
* OpenAddress database cleanup is aggressive
* Add backups for base data download and flood mapping (archived code)

application
* shut down shiny server on idle fix

printing
* Title page references featurefile (does not exist with aoi calls)