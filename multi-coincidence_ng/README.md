# Multi-coincidence neutron and gamma-ray detection

## About

This tally allows event-by-event "list-mode" collection of events with arbitrarily-set minimum number of required interaction coincidences (as inputted and then evaluated by a [Counter]) each satisfying a set minimum energy deposition threshold; 
these requirements can be specified individually for both neutron interactions and gamma-ray interactions. An arbitrarily long list of cells/regions treated within the coincidence logic can be provided.

This tally can be thought as something similar to a combination of the [T-Deposit2] (but with the number of detectors involved set arbitrarily rather than limited to 2) and [T-Product] tallies, producing something similar to a 
"dump" file but with "history counters" `chmax(1)` and `chmax(2)` allowed and enforced.

This tally was developed with a "camera-style" detector in mind that requires, for particle imaging, neutrons to undergo two interactions and gamma rays three interactions, each in a different region (scintillator bar) of the geometry (detector array). 
It outputs, for each coincidence event satisfying scoring criteria, the region number, energy deposition, interaction coordinates, and time of each interaction, along with additional identifying information and (optionally) extra information about the nature of the interactions. 
In this way, it minimally emulates the output of a real detector system (interaction energy deposition and location), though without the associated experimental systematic uncertainties, necessary for particle imaging 
but also allows a much more detailed look at what exact interactions contribute to the energy deposition in each interaction only available in simulation.

## PHITS Input

### [Counter] section

You should build a [Counter] section whose first counter increments by 1 when a neutron collision occurs in any of the "detector" cells and whose second counter similarly increments when a gamma ray collides in such cells.  An example is shown below:

```
$ Set counter to increment when:
$ neutron or photon undergo collision in a detector region
[Counter]  
counter = 1
part = neutron 
reg in     out    coll
100 10000  10000  0
200 0      0      1
210 0      0      1
220 0      0      1
230 0      0      1
240 0      0      1
900 0      0      0
counter = 2
part = photon 
reg in     out    coll
100 10000  10000  0
200 0      0      1
210 0      0      1
220 0      0      1
230 0      0      1
240 0      0      1
900 0      0      0
```

In this example, cells 200, 210, 220, 230, and 240 are "detector region" cells included in the coincidence logic.

### [T-Userdefined] section

In the [T-Userdefined] tally section in a PHITS input, values should be set as follows:

- **`nudtvar`** should be either (a) the number of cells `N` included in the coincidence logic or (b) this number of cells included `N` plus 6 (in which case, additional settings can be provided, detailed below)

In either case, the following first `N` entries for `udtvar`, `udtvar(1)` through `udtvar(N)`, should each be provided one cell number of a "detector" region involved in the coincidence logic.

Then, if `nudtvar` was set to `N` + 6, the following can also be specified:
- `udtvar(N+1)` / `udtvar(nudtvar-5)` = number of regions with above threshold neutron energy deposition required to be scored, "required `chmax(1)`"
- `udtvar(N+2)` / `udtvar(nudtvar-4)` = number of regions with above threshold gamma-ray energy deposition required to be scored, "required `chmax(2)`"
- `udtvar(N+3)` / `udtvar(nudtvar-3)` = minimum energy deposit threshold (in MeV) for neutrons to be recorded; if not exceeded, event is considered "not detected/detectable"
- `udtvar(N+4)` / `udtvar(nudtvar-2)` = minimum energy deposit threshold (in MeV) for gamma rays to be recorded; if not exceeded, event is considered "not detected/detectable"
- `udtvar(N+5)` / `udtvar(nudtvar-1)` = toggle for writing of extra informational reaction lines (`1`=yes, `0`=no), more info below
- `udtvar(N+6)` / `udtvar(nudtvar)` = `1` if wanting to specify the five above parameters as (b) earlier described; NOT `1` if using option (a) where only cell numbers are provided.

To be clear, the FINAL `udtvar(nudtvar)` specifies whether the five before it,  `(nudtvar-1)` to `(nudtvar-5)` mean to control the max counter values, energy deposition thresholds, and writing of extra reaction info or are just cell numbers (including this final variable).

If this final udtvar is not equal to `1`, then all `udtvar` values are interpreted to be cell numbers, the required `chmax(1)` and `chmax(2)` are assumed to be `1` (no multi-hit coincidence required), no energy threshold is enforced (E_deposit_min = 0), and extra reaction data lines are written.
*Note: If (b) and one of your cell numbers is `1`, it **must not** be provided to `udtvar(N)`; it instead should be provided to one of the other `udtvar(<N)` entries.*

If equal to one, then the preceeding five are for the desired values of the specified counters and control for whether the reaction lines, which describe secondary reaction products and incident particle energy, are written in addition to the standard energy deposition lines.


Following the same example as above, the [T-Userdefined] section could be written as shown below. 
In this example, neutron events require 2 different cells to experience energy depositions of at least 0.1 MeV, 
and gamma-ray events require at least 3 different cells each experiencing at least 0.05 MeV energy deposition in them.

```
[ T-Userdefined ]
nudtvar = 11 # should be either the number of cells or this number plus 6; here, 5+6=11
udtvar(1) = 200
udtvar(2) = 210
udtvar(3) = 220
udtvar(4) = 230
udtvar(5) = 240
# last six variables have different meaning 
udtvar( 6) = 2    # required final value of counter 1 / n colisions in different cells
udtvar( 7) = 3    # required final value of counter 2 / g colisions in different cells
# more accurately, at least one counter must be >0, these values are the total number
# of discrete regions with energy depositions exceeding the threshold values set below.
udtvar( 8) = 0.1  # minimum E deposit for neutron reactions to be recorded (MeV)
udtvar( 9) = 0.05 # minimum E deposit for gamma reactions to be recorded (MeV)
udtvar(10) = 1    # 1 write reaction lines too / 0 write only energy deposition lines
udtvar(11) = 1    # =1 specifies if above five are specified and used, otherwise means they don't exist
file = usrdef.out
```

## Output

The output filename is specified with the [T-Userdefined] `file` parameter. 
This file consists of event-by-event data for each set of coincident interactions within a history. 

After the header block displaying column header information, each coincident event registered will begin with either 
`ne` for neutron event or `ge` for gamma event, followed by the region number, energy deposition, spatial xyz coordinates, and time of each coincident reaction.

If additional reaction information writing is selected/enabled, immediately following each `ne`/`ge` line is at least one `rx` reaction line including information for each reaction. 
This line's reaction information, in order, is: `ncol`, `mathz`, `mathn`, `jcoll`, `kcoll`, `nclsts`, `nocas`, `no`, `name`, cell number, x coordinate, y coordinate, and z coordinate. 
It is then followed by an `In` line with the kf-code, incident energy, and incident particle weight of the particle initiating the reaction.
This is then followed by an `Out` line with the kf-code, energy, and particle weight of the produced particle(s), delimited by commas `,`.

These `rx` blocks allow more accurate determination of **how** energy was deposited (recoil nucleus species, whether multiple collisions occured in a single detector region, etc.).

Example output, with extra reaction info writing enabled, is shown below:

```
!Required final counter 1 value =       2 ; Required final counter 2 value =       3
!       #iomp    #batch  #history       #no     #name        #reg  EdepA(MeV)      xA(cm)      yA(cm)      zA(cm)      tA(ns)        #reg  EdepB(MeV)      xB(cm)      yB(cm)      zB(cm)      tB(ns)        #reg  EdepC(MeV)      xC(cm)      yC(cm)      zC(cm)      tC(ns)
!ncol   Z   N jcl kcl nclsts 
!In/Out kf-code     E(MeV)      weight 

ne          0         1     13202         2         7 ;       210  3.7104E-01  5.8454E+00  1.7331E+01  3.5470E-01  9.5980E+00 ,       220  1.7929E-01  1.6990E+00  1.8213E+01  3.8303E-01  1.3760E+01 ,
rx 14   6   6 10  3   2     13202         2         7         210              5.8454E+00  1.7331E+01  3.5470E-01
  In      2112  8.9615E-01  1.0000E+00
 Out      2112  7.9686E-01  1.0000E+00 ,   6000012  9.9337E-02  1.0000E+00
rx 14   1   0 10  3   2     13202         2         8         210              4.6594E+00  1.7611E+01 -2.4318E-01
  In      2112  7.9686E-01  1.0000E+00
 Out      2112  5.7381E-01  1.0000E+00 ,      2212  2.2325E-01  1.0000E+00
rx 14   1   0 10  3   2     13202         2         9         210              3.6476E+00  1.7987E+01 -1.4757E-01
  In      2112  5.7381E-01  1.0000E+00
 Out      2112  5.2551E-01  1.0000E+00 ,      2212  4.8448E-02  1.0000E+00
rx 14   6   6 10  3   2     13202         2        10         220              1.6990E+00  1.8213E+01  3.8303E-01
  In      2112  5.2551E-01  1.0000E+00
 Out      2112  5.2112E-01  1.0000E+00 ,   6000012  4.3940E-03  1.0000E+00
rx 14   1   0 10  3   2     13202         2        11         220              1.0280E+00  1.8433E+01  3.9513E-01
  In      2112  5.2112E-01  1.0000E+00
 Out      2112  3.4638E-01  1.0000E+00 ,      2212  1.7490E-01  1.0000E+00

ne          0         1     40165         7         3 ;       200  5.3470E-01  5.6593E+00  1.6461E+01 -1.5927E-01  1.0114E+01 ,       220  2.6337E-01  6.4522E+00  1.8655E+01 -1.7594E-01  1.1908E+01 ,
rx 14  14  15 10  3   2     40165         7         3         200              5.6593E+00  1.6461E+01 -1.5927E-01
  In      2112  1.4185E+00  1.0000E+00
 Out      2112  1.3952E+00  1.0000E+00 ,  14000029  2.3385E-02  1.0000E+00
rx 14   1   0 10  3   2     40165         7         4         200              5.6509E+00  1.6497E+01 -1.7044E-01
  In      2112  1.3952E+00  1.0000E+00
 Out      2112  8.8468E-01  1.0000E+00 ,      2212  5.1131E-01  1.0000E+00
rx 14   1   0 10  3   2     40165         7         5         220              6.4522E+00  1.8655E+01 -1.7594E-01
  In      2112  8.8468E-01  1.0000E+00
 Out      2112  6.2169E-01  1.0000E+00 ,      2212  2.6337E-01  1.0000E+00

ge          0         1    117867         6         1 ;       210  1.8221E+00  2.8781E+00  1.7140E+01 -1.2670E-01  3.8210E+00 ,       220  3.2727E+00  3.5527E+00  1.8200E+01 -4.7162E-02  3.8653E+00 ,       230  7.7301E-01  4.2937E+00  1.9600E+01 -2.6579E-01  3.9259E+00 ,
rx 14   0   0 13  0   3    117867         6         1         210              2.8781E+00  1.7140E+01 -1.2670E-01
  In        22  6.7861E+00  1.0000E+00
 Out        22  7.7937E-01  1.0000E+00 ,        11  6.0067E+00  1.0000E+00 ,        11  1.4000E-05  1.0000E+00
```


