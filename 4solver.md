
# 4 Solver Test

Below we will implement each of the 4 solvers unused on 4 different datasets. After separating out the best results of each solver, they will be merged together. Then we will pull the data from each solver on each of these merged solutions, to compare exactly how each result performed from each solver.

We will use the MEF2C gene as a base.

## Initialization

The first thing we do is load the sample data:



```R
# Load the assay matrix (mtx.sub) and transform it 3 ways
load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))

mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
mtx.log2 <- log2(mtx.tmp)

mtx.asinh <- asinh(mtx.sub)

suppressMessages(library(limma))
mtx.voom <- voom(mtx.sub)$E
```

## TReNA Data Analysis and Formatting

Now the function we use below does multiple things with our data. 
- First it runs the above code to create 4 distributions of data
- Then it runs TReNA on each distribution, using 4 different solvers
    - Spearman
    - Ridge
    - Lassopv
    - Ensemble (consisting of lasso, spearman, randomForest, ridge, lassopv, sqrtlasso, & pearson)
- That data is then organized in a single table (only consisting of the top 10 values of each test)
- The gene correlation data is added and the data is formatted properly, then returned!


```R
# Source the desired function and generate the table
suppressMessages(source("evalSolvers.R"))
tbl.all <- assess_ampAD154AllSolversAndDistributions()
```

    [1] "--- Testing Spearman"
    [1] "--- Testing Ridge"
    [1] "--- Testing Lassopv"
    [1] "--- Testing Ensemble"


Great! Now the data is stored in a convenient table:


```R
tbl.all
```


<table>
<thead><tr><th></th><th scope=col>gene</th><th scope=col>spearman.as.is</th><th scope=col>spearman.log2</th><th scope=col>spearman.asinh</th><th scope=col>spearman.voom</th><th scope=col>ridge.as.is</th><th scope=col>ridge.log2</th><th scope=col>ridge.asinh</th><th scope=col>ridge.voom</th><th scope=col>lpv.as.is</th><th scope=col>lpv.log2</th><th scope=col>lpv.asinh</th><th scope=col>lpv.voom</th><th scope=col>ensemble.as.is</th><th scope=col>ensemble.log2</th><th scope=col>ensemble.asinh</th><th scope=col>ensemble.voom</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>HLF          </td><td> 0.941550149 </td><td> 0.941550149 </td><td> 0.941550149 </td><td> 0.941697603 </td><td>  0.10740617 </td><td> 0.0421274219</td><td> 0.04472347  </td><td> 0.045332433 </td><td>52.742233460 </td><td>48.921247033 </td><td>48.38330171  </td><td>48.973631221 </td><td> 6.760785    </td><td> 6.817859    </td><td> 6.854995    </td><td> 6.967867    </td><td> 0.92323940  </td></tr>
	<tr><th scope=row>2</th><td>STAT4        </td><td> 0.920573708 </td><td> 0.920573708 </td><td> 0.920573708 </td><td> 0.922521887 </td><td>  0.59138702 </td><td> 0.0356272098</td><td> 0.03628105  </td><td> 0.038130485 </td><td>50.456971173 </td><td>54.221791193 </td><td>54.05992920  </td><td>54.624312937 </td><td>22.270523    </td><td>22.940225    </td><td>22.979982    </td><td>23.164242    </td><td> 0.90475987  </td></tr>
	<tr><th scope=row>3</th><td>SATB1        </td><td> 0.890366873 </td><td> 0.890366873 </td><td> 0.890366873 </td><td> 0.893133310 </td><td>  0.09068568 </td><td> 0.0495463228</td><td> 0.04830726  </td><td> 0.047961667 </td><td> 0.144372753 </td><td>42.768143862 </td><td> 0.01903403  </td><td> 0.002191662 </td><td> 2.662756    </td><td> 2.712079    </td><td> 2.694855    </td><td> 2.469122    </td><td> 0.86282704  </td></tr>
	<tr><th scope=row>4</th><td>SATB2        </td><td> 0.876561497 </td><td> 0.876561497 </td><td> 0.876561497 </td><td> 0.875618127 </td><td>  0.20449561 </td><td> 0.0545907748</td><td> 0.05485779  </td><td> 0.057232402 </td><td> 3.236844208 </td><td>48.154646517 </td><td>48.89733145  </td><td>49.424139376 </td><td> 2.348498    </td><td> 2.295249    </td><td> 2.407834    </td><td> 2.393848    </td><td> 0.84232989  </td></tr>
	<tr><th scope=row>5</th><td>ATF2         </td><td> 0.853196749 </td><td> 0.853196749 </td><td> 0.853196749 </td><td> 0.872160220 </td><td>  0.22520646 </td><td> 0.0761740549</td><td> 0.07325669  </td><td> 0.075470958 </td><td> 0.137274983 </td><td>21.089153375 </td><td> 0.29125507  </td><td> 0.135304746 </td><td> 2.008705    </td><td> 2.047409    </td><td> 2.019311    </td><td> 2.054066    </td><td> 0.82786857  </td></tr>
	<tr><th scope=row>6</th><td>FOXP2        </td><td> 0.831408752 </td><td> 0.831408752 </td><td> 0.831408752 </td><td> 0.836584162 </td><td>  1.20839645 </td><td> 0.0485015117</td><td> 0.05072036  </td><td> 0.053131039 </td><td>27.930161047 </td><td>30.743897428 </td><td>29.16804671  </td><td>29.826784860 </td><td>15.614550    </td><td>16.077252    </td><td>16.120648    </td><td>16.252081    </td><td> 0.81961721  </td></tr>
	<tr><th scope=row>8</th><td>TSHZ2        </td><td> 0.814784995 </td><td> 0.814784995 </td><td> 0.814784995 </td><td> 0.819328141 </td><td>  0.72360223 </td><td> 0.0473423421</td><td> 0.04597392  </td><td> 0.044437976 </td><td>12.605970976 </td><td>37.899825409 </td><td>37.06080200  </td><td>37.505220279 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.77043615  </td></tr>
	<tr><th scope=row>7</th><td>DRGX         </td><td> 0.816968343 </td><td> 0.816968343 </td><td> 0.816968343 </td><td> 0.819446551 </td><td>  3.31900405 </td><td> 0.0156392744</td><td> 0.02333311  </td><td> 0.028977456 </td><td> 0.043912086 </td><td> 0.141846140 </td><td> 0.06635882  </td><td> 0.016631006 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.73094500  </td></tr>
	<tr><th scope=row>10</th><td>SOX13        </td><td>-0.767238358 </td><td>-0.767238358 </td><td>-0.767238358 </td><td>-0.761502734 </td><td> -0.25073994 </td><td>-0.0349220623</td><td>-0.03368853  </td><td>-0.036125864 </td><td>12.408198928 </td><td>16.017307839 </td><td>13.70820824  </td><td>15.936759724 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.72731213  </td></tr>
	<tr><th scope=row>14</th><td>HDX          </td><td> 0.701936366 </td><td> 0.701936366 </td><td> 0.701936366 </td><td> 0.721930340 </td><td>  4.14010657 </td><td> 0.0474368599</td><td> 0.04498817  </td><td> 0.053274039 </td><td>27.446268097 </td><td> 1.404191761 </td><td> 0.52651363  </td><td> 0.591778019 </td><td>84.112576    </td><td>86.693371    </td><td>86.919009    </td><td>87.605358    </td><td> 0.72419056  </td></tr>
	<tr><th scope=row>13</th><td>LHX6         </td><td> 0.728602631 </td><td> 0.728602631 </td><td> 0.728602631 </td><td> 0.742710176 </td><td>  0.79829553 </td><td> 0.0423595791</td><td> 0.04224198  </td><td> 0.042939693 </td><td>20.055971685 </td><td>24.763932878 </td><td>26.04871481  </td><td>25.916705889 </td><td>17.284988    </td><td>17.982962    </td><td>17.935697    </td><td>18.287767    </td><td> 0.71461154  </td></tr>
	<tr><th scope=row>9</th><td>ESRRG        </td><td> 0.813791915 </td><td> 0.813791915 </td><td> 0.813791915 </td><td> 0.818304341 </td><td>  0.70051147 </td><td> 0.0395895408</td><td> 0.04189497  </td><td> 0.043570742 </td><td> 3.904288201 </td><td> 2.816872225 </td><td> 2.82494002  </td><td> 1.607153573 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.70903181  </td></tr>
	<tr><th scope=row>12</th><td>TSHZ3        </td><td> 0.745134229 </td><td> 0.745134229 </td><td> 0.745134229 </td><td> 0.756719976 </td><td>  0.66204547 </td><td> 0.0726322744</td><td> 0.07453427  </td><td> 0.074190934 </td><td> 0.036665968 </td><td>20.851479709 </td><td>27.51650187  </td><td>29.734229294 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.70562897  </td></tr>
	<tr><th scope=row>16</th><td>FOXO3        </td><td> 0.669744489 </td><td> 0.669744489 </td><td> 0.669744489 </td><td> 0.776705572 </td><td>  0.26750560 </td><td> 0.0583238597</td><td> 0.05719985  </td><td> 0.054122714 </td><td> 0.059053565 </td><td> 0.110589663 </td><td> 0.16735094  </td><td> 0.288591705 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.64598153  </td></tr>
	<tr><th scope=row>15</th><td>NFE2L2       </td><td>-0.672528240 </td><td>-0.672528240 </td><td>-0.672528240 </td><td>-0.708746616 </td><td> -0.11601251 </td><td>-0.0267063552</td><td>-0.02391370  </td><td>-0.027447909 </td><td> 9.438418560 </td><td> 4.016293261 </td><td> 2.44841020  </td><td> 6.954817011 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.63838196  </td></tr>
	<tr><th scope=row>17</th><td>POU3F3       </td><td>-0.668608423 </td><td>-0.668608423 </td><td>-0.668608423 </td><td>-0.687904224 </td><td> -0.08553820 </td><td>-0.0398283814</td><td>-0.03996137  </td><td>-0.057652276 </td><td> 0.045183852 </td><td> 0.188709445 </td><td> 0.14182217  </td><td> 1.211854479 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.62965384  </td></tr>
	<tr><th scope=row>18</th><td>SOX12        </td><td>-0.667379640 </td><td>-0.667379640 </td><td>-0.667379640 </td><td>-0.639926396 </td><td> -0.34395997 </td><td>-0.0564014679</td><td>-0.05354446  </td><td>-0.060381157 </td><td>15.918749985 </td><td>18.646102685 </td><td>20.19427473  </td><td>17.964569412 </td><td> 1.236455    </td><td> 1.243567    </td><td> 1.255499    </td><td> 1.240401    </td><td>-0.62404571  </td></tr>
	<tr><th scope=row>11</th><td>FOXO4        </td><td>-0.753055969 </td><td>-0.753055969 </td><td>-0.753055969 </td><td>-0.773118084 </td><td> -0.01320554 </td><td>-0.0270975633</td><td>-0.02807569  </td><td>-0.035534301 </td><td> 0.125579014 </td><td>28.913023119 </td><td>27.17059677  </td><td>30.623290540 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.55645498  </td></tr>
	<tr><th scope=row>21</th><td>FOXK2        </td><td> 0.551600923 </td><td> 0.551600923 </td><td> 0.551600923 </td><td> 0.590891256 </td><td>  0.20327471 </td><td> 0.0731278625</td><td> 0.07325921  </td><td> 0.072778810 </td><td> 0.244015522 </td><td> 0.391873801 </td><td> 0.39598950  </td><td> 0.269501342 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.53232915  </td></tr>
	<tr><th scope=row>20</th><td>FOXP1        </td><td> 0.589628961 </td><td> 0.589628961 </td><td> 0.589628961 </td><td> 0.615828288 </td><td>  0.27595213 </td><td> 0.0754438536</td><td> 0.08134690  </td><td> 0.075309257 </td><td> 2.365100032 </td><td> 0.894125005 </td><td> 1.27649028  </td><td> 0.232164762 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.51591782  </td></tr>
	<tr><th scope=row>22</th><td>TAF1         </td><td>-0.518577941 </td><td>-0.518577941 </td><td>-0.518577941 </td><td>-0.346293083 </td><td> -0.43906196 </td><td>-0.0612173561</td><td>-0.05580276  </td><td>-0.035123372 </td><td> 0.429184476 </td><td> 0.263192726 </td><td> 0.15618020  </td><td> 0.085302124 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.47800023  </td></tr>
	<tr><th scope=row>19</th><td>STAT1        </td><td> 0.656799819 </td><td> 0.656799819 </td><td> 0.656799819 </td><td> 0.674345163 </td><td>  0.09303440 </td><td> 0.0435315637</td><td> 0.04073636  </td><td> 0.038480909 </td><td>10.942170686 </td><td>14.926963003 </td><td>15.51621561  </td><td>10.775635153 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.44791101  </td></tr>
	<tr><th scope=row>24</th><td>ADNP2        </td><td> 0.421795090 </td><td> 0.421795090 </td><td> 0.421795090 </td><td> 0.484999772 </td><td>  0.97393424 </td><td> 0.0670899314</td><td> 0.06529808  </td><td> 0.061587778 </td><td> 3.439154051 </td><td> 0.003718927 </td><td> 0.08054365  </td><td> 0.132897214 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.43089283  </td></tr>
	<tr><th scope=row>27</th><td>NR3C1        </td><td> 0.351312103 </td><td> 0.351312103 </td><td> 0.351312103 </td><td> 0.413744886 </td><td>  0.26035196 </td><td> 0.0556031269</td><td> 0.05495463  </td><td> 0.060177477 </td><td> 0.080191130 </td><td> 1.038078556 </td><td> 0.01149814  </td><td> 0.035837630 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.40598438  </td></tr>
	<tr><th scope=row>23</th><td>STAT5B       </td><td>-0.422379320 </td><td>-0.422379320 </td><td>-0.422379320 </td><td>-0.432139766 </td><td> -0.17967771 </td><td>-0.0930294904</td><td>-0.09419464  </td><td>-0.138263144 </td><td> 0.388318520 </td><td> 3.436162729 </td><td> 5.14780111  </td><td> 9.666864945 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.37452622  </td></tr>
	<tr><th scope=row>25</th><td>NR5A2        </td><td>-0.404429112 </td><td>-0.404429112 </td><td>-0.404429112 </td><td>-0.278828434 </td><td> -6.38009764 </td><td> 0.0005382409</td><td>-0.04530581  </td><td>-0.020322057 </td><td> 0.743878327 </td><td> 0.158231009 </td><td> 0.78869852  </td><td> 0.602078521 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.36411016  </td></tr>
	<tr><th scope=row>26</th><td>FOXD4L3      </td><td> 0.396624885 </td><td> 0.396624885 </td><td> 0.396624885 </td><td> 0.437130300 </td><td>  7.15764709 </td><td> 0.0072629714</td><td> 0.03595524  </td><td> 0.032994383 </td><td> 0.274378117 </td><td> 0.079173770 </td><td> 0.17540019  </td><td> 0.061391763 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.29528836  </td></tr>
	<tr><th scope=row>28</th><td>RAX2         </td><td>-0.345922104 </td><td>-0.345922104 </td><td>-0.345922104 </td><td>-0.185103018 </td><td>-13.39956056 </td><td>-0.0035605472</td><td>-0.08226796  </td><td>-0.047156084 </td><td> 0.273495549 </td><td> 0.058783176 </td><td> 5.88962964  </td><td> 5.310608493 </td><td>26.091555    </td><td>25.890837    </td><td>26.462607    </td><td>25.662322    </td><td>-0.28106544  </td></tr>
	<tr><th scope=row>29</th><td>POU5F2       </td><td>-0.222980251 </td><td>-0.222980251 </td><td>-0.222980251 </td><td>-0.114969676 </td><td>-13.48946528 </td><td>-0.0055976068</td><td>-0.01666904  </td><td> 0.006172083 </td><td> 0.978019558 </td><td> 0.094400101 </td><td> 0.08626123  </td><td> 0.212442690 </td><td>13.743087    </td><td>13.727787    </td><td>13.773426    </td><td>13.768085    </td><td>-0.27661563  </td></tr>
	<tr><th scope=row>31</th><td>FOXE3        </td><td>-0.141698283 </td><td>-0.141698283 </td><td>-0.141698283 </td><td> 0.014661335 </td><td>-18.27342871 </td><td>-0.0002076319</td><td>-0.06017492  </td><td>-0.012819682 </td><td> 1.705123299 </td><td> 0.015505164 </td><td> 1.46269473  </td><td> 0.065927361 </td><td>18.398443    </td><td>18.393907    </td><td>18.324356    </td><td>18.405406    </td><td>-0.15567610  </td></tr>
	<tr><th scope=row>35</th><td>FOXD4        </td><td>-0.042989808 </td><td>-0.042989808 </td><td>-0.042989808 </td><td> 0.009325066 </td><td> -7.22056969 </td><td>-0.0049860013</td><td>-0.02327859  </td><td>-0.013548884 </td><td> 5.382415210 </td><td> 1.082873748 </td><td> 2.50618718  </td><td> 1.216426098 </td><td>10.955887    </td><td>10.809892    </td><td>11.010188    </td><td>10.677929    </td><td>-0.11598769  </td></tr>
	<tr><th scope=row>30</th><td>OTX2         </td><td>-0.146096961 </td><td>-0.146096961 </td><td>-0.146096961 </td><td> 0.111551984 </td><td> -0.60068988 </td><td>-0.0028640186</td><td>-0.07130008  </td><td>-0.057415723 </td><td> 0.074321209 </td><td> 0.040210069 </td><td> 3.68163888  </td><td> 2.487164341 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.11114757  </td></tr>
	<tr><th scope=row>34</th><td>UNCX         </td><td>-0.095060793 </td><td>-0.095060793 </td><td>-0.095060793 </td><td> 0.206505623 </td><td> -0.50063456 </td><td>-0.0064316424</td><td>-0.06454465  </td><td>-0.048802587 </td><td> 0.001569643 </td><td> 1.703867927 </td><td> 0.03213779  </td><td> 0.012959038 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.10945671  </td></tr>
	<tr><th scope=row>32</th><td>FOXP3        </td><td>-0.130418269 </td><td>-0.130418269 </td><td>-0.130418269 </td><td>-0.016774841 </td><td>  5.63388607 </td><td> 0.0101754298</td><td> 0.01456269  </td><td> 0.022628827 </td><td> 0.781730811 </td><td> 1.963567382 </td><td> 0.88734673  </td><td> 0.846382666 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.10724617  </td></tr>
	<tr><th scope=row>33</th><td>ELF3         </td><td>-0.111150147 </td><td>-0.111150147 </td><td>-0.111150147 </td><td> 0.020903552 </td><td> -5.78714019 </td><td>-0.0020871648</td><td>-0.03818554  </td><td>-0.002721941 </td><td> 0.217839610 </td><td> 0.080522787 </td><td> 0.16783590  </td><td> 0.160927583 </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td> 0.000000    </td><td>-0.10576117  </td></tr>
	<tr><th scope=row>37</th><td>FOXD4L1      </td><td>-0.007069132 </td><td>-0.007069132 </td><td>-0.007069132 </td><td> 0.014811582 </td><td>-10.77880221 </td><td>-0.0072751002</td><td>-0.04523622  </td><td>-0.036691629 </td><td> 4.644059166 </td><td> 1.321882550 </td><td> 2.89433216  </td><td> 2.771612886 </td><td>17.852311    </td><td>17.681981    </td><td>17.976983    </td><td>17.484736    </td><td>-0.10509003  </td></tr>
	<tr><th scope=row>36</th><td>IRF4         </td><td>-0.014939560 </td><td>-0.014939560 </td><td>-0.014939560 </td><td> 0.109388767 </td><td>  9.70058260 </td><td> 0.0031461546</td><td> 0.05225426  </td><td> 0.034519887 </td><td> 0.294906912 </td><td> 0.351799360 </td><td> 0.96531799  </td><td> 0.251812993 </td><td> 9.960390    </td><td> 9.952887    </td><td> 9.956741    </td><td> 9.971431    </td><td> 0.01753765  </td></tr>
</tbody>
</table>



Using the commands below, we can quickly plot the results from the table and compare each of the 4 solvers. The results are especially interesting when you note that these plots were generated from the same datasets.


```R
par(family = "sans")
pairs(tbl.all[2:18],
     labels = names(tbl.all)[2:18])
```


![png](output_8_0.png)


## Individual Solver Performance and Correlation

### Spearman

These are a ton of graphs to look at all at once, and each one is packed with information. Let's take a look at how each solver performed over the different datasets, starting with spearman.


```R
head(tbl.all[,c(1,2,3,4,5,18)],10)
```


<table>
<thead><tr><th></th><th scope=col>gene</th><th scope=col>spearman.as.is</th><th scope=col>spearman.log2</th><th scope=col>spearman.asinh</th><th scope=col>spearman.voom</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>HLF       </td><td> 0.9415501</td><td> 0.9415501</td><td> 0.9415501</td><td> 0.9416976</td><td> 0.9232394</td></tr>
	<tr><th scope=row>2</th><td>STAT4     </td><td> 0.9205737</td><td> 0.9205737</td><td> 0.9205737</td><td> 0.9225219</td><td> 0.9047599</td></tr>
	<tr><th scope=row>3</th><td>SATB1     </td><td> 0.8903669</td><td> 0.8903669</td><td> 0.8903669</td><td> 0.8931333</td><td> 0.8628270</td></tr>
	<tr><th scope=row>4</th><td>SATB2     </td><td> 0.8765615</td><td> 0.8765615</td><td> 0.8765615</td><td> 0.8756181</td><td> 0.8423299</td></tr>
	<tr><th scope=row>5</th><td>ATF2      </td><td> 0.8531967</td><td> 0.8531967</td><td> 0.8531967</td><td> 0.8721602</td><td> 0.8278686</td></tr>
	<tr><th scope=row>6</th><td>FOXP2     </td><td> 0.8314088</td><td> 0.8314088</td><td> 0.8314088</td><td> 0.8365842</td><td> 0.8196172</td></tr>
	<tr><th scope=row>8</th><td>TSHZ2     </td><td> 0.8147850</td><td> 0.8147850</td><td> 0.8147850</td><td> 0.8193281</td><td> 0.7704362</td></tr>
	<tr><th scope=row>7</th><td>DRGX      </td><td> 0.8169683</td><td> 0.8169683</td><td> 0.8169683</td><td> 0.8194466</td><td> 0.7309450</td></tr>
	<tr><th scope=row>10</th><td>SOX13     </td><td>-0.7672384</td><td>-0.7672384</td><td>-0.7672384</td><td>-0.7615027</td><td>-0.7273121</td></tr>
	<tr><th scope=row>14</th><td>HDX       </td><td> 0.7019364</td><td> 0.7019364</td><td> 0.7019364</td><td> 0.7219303</td><td> 0.7241906</td></tr>
</tbody>
</table>



Okay, this table is a bit easier to digest. From a cursory glance, it looks as though spearman returned nearly the exact same outputs for each dataset! A closer look proves that these results are very similar.


```R
cor(tbl.all[2:5])
```


<table>
<thead><tr><th></th><th scope=col>spearman.as.is</th><th scope=col>spearman.log2</th><th scope=col>spearman.asinh</th><th scope=col>spearman.voom</th></tr></thead>
<tbody>
	<tr><th scope=row>spearman.as.is</th><td>1.0000000</td><td>1.0000000</td><td>1.0000000</td><td>0.9911316</td></tr>
	<tr><th scope=row>spearman.log2</th><td>1.0000000</td><td>1.0000000</td><td>1.0000000</td><td>0.9911316</td></tr>
	<tr><th scope=row>spearman.asinh</th><td>1.0000000</td><td>1.0000000</td><td>1.0000000</td><td>0.9911316</td></tr>
	<tr><th scope=row>spearman.voom</th><td>0.9911316</td><td>0.9911316</td><td>0.9911316</td><td>1.0000000</td></tr>
</tbody>
</table>



In this table we ignored the specific genes and gene correlation to just compare data directly from the solver. The largest difference between any of these results is less than 1%, very close. Now we'll take a look at ridge's performance and similarity.

### Ridge


```R
head(tbl.all[,c(1,6,7,8,9,18)],10)
```


<table>
<thead><tr><th></th><th scope=col>gene</th><th scope=col>ridge.as.is</th><th scope=col>ridge.log2</th><th scope=col>ridge.asinh</th><th scope=col>ridge.voom</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>HLF        </td><td> 0.10740617</td><td> 0.04212742</td><td> 0.04472347</td><td> 0.04533243</td><td> 0.9232394 </td></tr>
	<tr><th scope=row>2</th><td>STAT4      </td><td> 0.59138702</td><td> 0.03562721</td><td> 0.03628105</td><td> 0.03813049</td><td> 0.9047599 </td></tr>
	<tr><th scope=row>3</th><td>SATB1      </td><td> 0.09068568</td><td> 0.04954632</td><td> 0.04830726</td><td> 0.04796167</td><td> 0.8628270 </td></tr>
	<tr><th scope=row>4</th><td>SATB2      </td><td> 0.20449561</td><td> 0.05459077</td><td> 0.05485779</td><td> 0.05723240</td><td> 0.8423299 </td></tr>
	<tr><th scope=row>5</th><td>ATF2       </td><td> 0.22520646</td><td> 0.07617405</td><td> 0.07325669</td><td> 0.07547096</td><td> 0.8278686 </td></tr>
	<tr><th scope=row>6</th><td>FOXP2      </td><td> 1.20839645</td><td> 0.04850151</td><td> 0.05072036</td><td> 0.05313104</td><td> 0.8196172 </td></tr>
	<tr><th scope=row>8</th><td>TSHZ2      </td><td> 0.72360223</td><td> 0.04734234</td><td> 0.04597392</td><td> 0.04443798</td><td> 0.7704362 </td></tr>
	<tr><th scope=row>7</th><td>DRGX       </td><td> 3.31900405</td><td> 0.01563927</td><td> 0.02333311</td><td> 0.02897746</td><td> 0.7309450 </td></tr>
	<tr><th scope=row>10</th><td>SOX13      </td><td>-0.25073994</td><td>-0.03492206</td><td>-0.03368853</td><td>-0.03612586</td><td>-0.7273121 </td></tr>
	<tr><th scope=row>14</th><td>HDX        </td><td> 4.14010657</td><td> 0.04743686</td><td> 0.04498817</td><td> 0.05327404</td><td> 0.7241906 </td></tr>
</tbody>
</table>



These data definitely vary more than the spearman data. Once again, we'll look at exactly how much they vary.


```R
cor(tbl.all[6:9])
```


<table>
<thead><tr><th></th><th scope=col>ridge.as.is</th><th scope=col>ridge.log2</th><th scope=col>ridge.asinh</th><th scope=col>ridge.voom</th></tr></thead>
<tbody>
	<tr><th scope=row>ridge.as.is</th><td>1.0000000</td><td>0.2092563</td><td>0.5245586</td><td>0.3389476</td></tr>
	<tr><th scope=row>ridge.log2</th><td>0.2092563</td><td>1.0000000</td><td>0.8779865</td><td>0.9340204</td></tr>
	<tr><th scope=row>ridge.asinh</th><td>0.5245586</td><td>0.8779865</td><td>1.0000000</td><td>0.9533540</td></tr>
	<tr><th scope=row>ridge.voom</th><td>0.3389476</td><td>0.9340204</td><td>0.9533540</td><td>1.0000000</td></tr>
</tbody>
</table>



There are some recognizeable patterns in this table. Log2, voom, and asinh (but especially the first two) are pretty similar, while the as-is dataset resulted in very different outputs. Onto lassopv:

### Lassopv


```R
head(tbl.all[,c(1,10,11,12,13,18)],10)
cor(tbl.all[10:13])
```


<table>
<thead><tr><th></th><th scope=col>gene</th><th scope=col>lpv.as.is</th><th scope=col>lpv.log2</th><th scope=col>lpv.asinh</th><th scope=col>lpv.voom</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>HLF         </td><td>52.74223346 </td><td>48.9212470  </td><td>48.38330171 </td><td>48.973631221</td><td> 0.9232394  </td></tr>
	<tr><th scope=row>2</th><td>STAT4       </td><td>50.45697117 </td><td>54.2217912  </td><td>54.05992920 </td><td>54.624312937</td><td> 0.9047599  </td></tr>
	<tr><th scope=row>3</th><td>SATB1       </td><td> 0.14437275 </td><td>42.7681439  </td><td> 0.01903403 </td><td> 0.002191662</td><td> 0.8628270  </td></tr>
	<tr><th scope=row>4</th><td>SATB2       </td><td> 3.23684421 </td><td>48.1546465  </td><td>48.89733145 </td><td>49.424139376</td><td> 0.8423299  </td></tr>
	<tr><th scope=row>5</th><td>ATF2        </td><td> 0.13727498 </td><td>21.0891534  </td><td> 0.29125507 </td><td> 0.135304746</td><td> 0.8278686  </td></tr>
	<tr><th scope=row>6</th><td>FOXP2       </td><td>27.93016105 </td><td>30.7438974  </td><td>29.16804671 </td><td>29.826784860</td><td> 0.8196172  </td></tr>
	<tr><th scope=row>8</th><td>TSHZ2       </td><td>12.60597098 </td><td>37.8998254  </td><td>37.06080200 </td><td>37.505220279</td><td> 0.7704362  </td></tr>
	<tr><th scope=row>7</th><td>DRGX        </td><td> 0.04391209 </td><td> 0.1418461  </td><td> 0.06635882 </td><td> 0.016631006</td><td> 0.7309450  </td></tr>
	<tr><th scope=row>10</th><td>SOX13       </td><td>12.40819893 </td><td>16.0173078  </td><td>13.70820824 </td><td>15.936759724</td><td>-0.7273121  </td></tr>
	<tr><th scope=row>14</th><td>HDX         </td><td>27.44626810 </td><td> 1.4041918  </td><td> 0.52651363 </td><td> 0.591778019</td><td> 0.7241906  </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>lpv.as.is</th><th scope=col>lpv.log2</th><th scope=col>lpv.asinh</th><th scope=col>lpv.voom</th></tr></thead>
<tbody>
	<tr><th scope=row>lpv.as.is</th><td>1.0000000</td><td>0.6274871</td><td>0.6970330</td><td>0.6842001</td></tr>
	<tr><th scope=row>lpv.log2</th><td>0.6274871</td><td>1.0000000</td><td>0.8799655</td><td>0.8769244</td></tr>
	<tr><th scope=row>lpv.asinh</th><td>0.6970330</td><td>0.8799655</td><td>1.0000000</td><td>0.9947926</td></tr>
	<tr><th scope=row>lpv.voom</th><td>0.6842001</td><td>0.8769244</td><td>0.9947926</td><td>1.0000000</td></tr>
</tbody>
</table>



It's difficult to tell from the table exactly what the correlations between each result are. The data does seem to vary wildly, from values less than 1 to near 50.
Although the data seemed to differ greatly, the correlations are actually better than much of ridge's as-is correlations. The final analysis will be of ensemble's datasets.

### Ensemble


```R
head(tbl.all[,c(1,14:18)],10)
cor(tbl.all[14:17])
```


<table>
<thead><tr><th></th><th scope=col>gene</th><th scope=col>ensemble.as.is</th><th scope=col>ensemble.log2</th><th scope=col>ensemble.asinh</th><th scope=col>ensemble.voom</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>HLF       </td><td> 6.760785 </td><td> 6.817859 </td><td> 6.854995 </td><td> 6.967867 </td><td> 0.9232394</td></tr>
	<tr><th scope=row>2</th><td>STAT4     </td><td>22.270523 </td><td>22.940225 </td><td>22.979982 </td><td>23.164242 </td><td> 0.9047599</td></tr>
	<tr><th scope=row>3</th><td>SATB1     </td><td> 2.662756 </td><td> 2.712079 </td><td> 2.694855 </td><td> 2.469122 </td><td> 0.8628270</td></tr>
	<tr><th scope=row>4</th><td>SATB2     </td><td> 2.348498 </td><td> 2.295249 </td><td> 2.407834 </td><td> 2.393848 </td><td> 0.8423299</td></tr>
	<tr><th scope=row>5</th><td>ATF2      </td><td> 2.008705 </td><td> 2.047409 </td><td> 2.019311 </td><td> 2.054066 </td><td> 0.8278686</td></tr>
	<tr><th scope=row>6</th><td>FOXP2     </td><td>15.614550 </td><td>16.077252 </td><td>16.120648 </td><td>16.252081 </td><td> 0.8196172</td></tr>
	<tr><th scope=row>8</th><td>TSHZ2     </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.7704362</td></tr>
	<tr><th scope=row>7</th><td>DRGX      </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.7309450</td></tr>
	<tr><th scope=row>10</th><td>SOX13     </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.000000 </td><td> 0.000000 </td><td>-0.7273121</td></tr>
	<tr><th scope=row>14</th><td>HDX       </td><td>84.112576 </td><td>86.693371 </td><td>86.919009 </td><td>87.605358 </td><td> 0.7241906</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>ensemble.as.is</th><th scope=col>ensemble.log2</th><th scope=col>ensemble.asinh</th><th scope=col>ensemble.voom</th></tr></thead>
<tbody>
	<tr><th scope=row>ensemble.as.is</th><td>1.0000000</td><td>0.9998912</td><td>0.9999406</td><td>0.9997643</td></tr>
	<tr><th scope=row>ensemble.log2</th><td>0.9998912</td><td>1.0000000</td><td>0.9999796</td><td>0.9999707</td></tr>
	<tr><th scope=row>ensemble.asinh</th><td>0.9999406</td><td>0.9999796</td><td>1.0000000</td><td>0.9999146</td></tr>
	<tr><th scope=row>ensemble.voom</th><td>0.9997643</td><td>0.9999707</td><td>0.9999146</td><td>1.0000000</td></tr>
</tbody>
</table>



This data seems very correlated. Note that 4 of the genes have NA data (which we filled in with 0's so they wouldn't break our calculations). Indeed, the data is very closely correlated. It's almost as closely correlated as spearman's data that we initially took a look at!

## Inter-Solver Analysis

Now that we've assessed how similar the data is between the datasets within each solver, it's time for us to analyze whether any solvers functioned better than others on specific datasets.

### as-is Dataset


```R
head(tbl.all[,c(1,2,6,10,14,18)],10)
cor(tbl.all[c(2,6,10,14,18)])
```


<table>
<thead><tr><th></th><th scope=col>gene</th><th scope=col>spearman.as.is</th><th scope=col>ridge.as.is</th><th scope=col>lpv.as.is</th><th scope=col>ensemble.as.is</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>HLF        </td><td> 0.9415501 </td><td> 0.10740617</td><td>52.74223346</td><td> 6.760785  </td><td> 0.9232394 </td></tr>
	<tr><th scope=row>2</th><td>STAT4      </td><td> 0.9205737 </td><td> 0.59138702</td><td>50.45697117</td><td>22.270523  </td><td> 0.9047599 </td></tr>
	<tr><th scope=row>3</th><td>SATB1      </td><td> 0.8903669 </td><td> 0.09068568</td><td> 0.14437275</td><td> 2.662756  </td><td> 0.8628270 </td></tr>
	<tr><th scope=row>4</th><td>SATB2      </td><td> 0.8765615 </td><td> 0.20449561</td><td> 3.23684421</td><td> 2.348498  </td><td> 0.8423299 </td></tr>
	<tr><th scope=row>5</th><td>ATF2       </td><td> 0.8531967 </td><td> 0.22520646</td><td> 0.13727498</td><td> 2.008705  </td><td> 0.8278686 </td></tr>
	<tr><th scope=row>6</th><td>FOXP2      </td><td> 0.8314088 </td><td> 1.20839645</td><td>27.93016105</td><td>15.614550  </td><td> 0.8196172 </td></tr>
	<tr><th scope=row>8</th><td>TSHZ2      </td><td> 0.8147850 </td><td> 0.72360223</td><td>12.60597098</td><td> 0.000000  </td><td> 0.7704362 </td></tr>
	<tr><th scope=row>7</th><td>DRGX       </td><td> 0.8169683 </td><td> 3.31900405</td><td> 0.04391209</td><td> 0.000000  </td><td> 0.7309450 </td></tr>
	<tr><th scope=row>10</th><td>SOX13      </td><td>-0.7672384 </td><td>-0.25073994</td><td>12.40819893</td><td> 0.000000  </td><td>-0.7273121 </td></tr>
	<tr><th scope=row>14</th><td>HDX        </td><td> 0.7019364 </td><td> 4.14010657</td><td>27.44626810</td><td>84.112576  </td><td> 0.7241906 </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>spearman.as.is</th><th scope=col>ridge.as.is</th><th scope=col>lpv.as.is</th><th scope=col>ensemble.as.is</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>spearman.as.is</th><td>1.0000000 </td><td> 0.3101389</td><td>0.3565977 </td><td> 0.1609107</td><td>0.9949736 </td></tr>
	<tr><th scope=row>ridge.as.is</th><td>0.3101389 </td><td> 1.0000000</td><td>0.1412356 </td><td>-0.1245656</td><td>0.3328061 </td></tr>
	<tr><th scope=row>lpv.as.is</th><td>0.3565977 </td><td> 0.1412356</td><td>1.0000000 </td><td> 0.4161600</td><td>0.3728764 </td></tr>
	<tr><th scope=row>ensemble.as.is</th><td>0.1609107 </td><td>-0.1245656</td><td>0.4161600 </td><td> 1.0000000</td><td>0.1776926 </td></tr>
	<tr><th scope=row>gene.cor</th><td>0.9949736 </td><td> 0.3328061</td><td>0.3728764 </td><td> 0.1776926</td><td>1.0000000 </td></tr>
</tbody>
</table>



Spearman by far outshines the rest of the solvers with the as-is data. Other than that, none of the other data seems to correlated very much.

### log2 Dataset


```R
head(tbl.all[,c(1,3,7,11,15,18)],10)
cor(tbl.all[c(3,7,11,15,18)])
```


<table>
<thead><tr><th></th><th scope=col>gene</th><th scope=col>spearman.log2</th><th scope=col>ridge.log2</th><th scope=col>lpv.log2</th><th scope=col>ensemble.log2</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>HLF        </td><td> 0.9415501 </td><td> 0.04212742</td><td>48.9212470 </td><td> 6.817859  </td><td> 0.9232394 </td></tr>
	<tr><th scope=row>2</th><td>STAT4      </td><td> 0.9205737 </td><td> 0.03562721</td><td>54.2217912 </td><td>22.940225  </td><td> 0.9047599 </td></tr>
	<tr><th scope=row>3</th><td>SATB1      </td><td> 0.8903669 </td><td> 0.04954632</td><td>42.7681439 </td><td> 2.712079  </td><td> 0.8628270 </td></tr>
	<tr><th scope=row>4</th><td>SATB2      </td><td> 0.8765615 </td><td> 0.05459077</td><td>48.1546465 </td><td> 2.295249  </td><td> 0.8423299 </td></tr>
	<tr><th scope=row>5</th><td>ATF2       </td><td> 0.8531967 </td><td> 0.07617405</td><td>21.0891534 </td><td> 2.047409  </td><td> 0.8278686 </td></tr>
	<tr><th scope=row>6</th><td>FOXP2      </td><td> 0.8314088 </td><td> 0.04850151</td><td>30.7438974 </td><td>16.077252  </td><td> 0.8196172 </td></tr>
	<tr><th scope=row>8</th><td>TSHZ2      </td><td> 0.8147850 </td><td> 0.04734234</td><td>37.8998254 </td><td> 0.000000  </td><td> 0.7704362 </td></tr>
	<tr><th scope=row>7</th><td>DRGX       </td><td> 0.8169683 </td><td> 0.01563927</td><td> 0.1418461 </td><td> 0.000000  </td><td> 0.7309450 </td></tr>
	<tr><th scope=row>10</th><td>SOX13      </td><td>-0.7672384 </td><td>-0.03492206</td><td>16.0173078 </td><td> 0.000000  </td><td>-0.7273121 </td></tr>
	<tr><th scope=row>14</th><td>HDX        </td><td> 0.7019364 </td><td> 0.04743686</td><td> 1.4041918 </td><td>86.693371  </td><td> 0.7241906 </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>spearman.log2</th><th scope=col>ridge.log2</th><th scope=col>lpv.log2</th><th scope=col>ensemble.log2</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>spearman.log2</th><td>1.0000000  </td><td>0.8562438  </td><td>0.459553521</td><td>0.165078751</td><td>0.9949736  </td></tr>
	<tr><th scope=row>ridge.log2</th><td>0.8562438  </td><td>1.0000000  </td><td>0.291736334</td><td>0.113554597</td><td>0.8613929  </td></tr>
	<tr><th scope=row>lpv.log2</th><td>0.4595535  </td><td>0.2917363  </td><td>1.000000000</td><td>0.001030541</td><td>0.4894849  </td></tr>
	<tr><th scope=row>ensemble.log2</th><td>0.1650788  </td><td>0.1135546  </td><td>0.001030541</td><td>1.000000000</td><td>0.1822591  </td></tr>
	<tr><th scope=row>gene.cor</th><td>0.9949736  </td><td>0.8613929  </td><td>0.489484901</td><td>0.182259127</td><td>1.0000000  </td></tr>
</tbody>
</table>



Once again, spearman has the best correlation, however ridge also excelled with the log2 data. 

### asinh Dataset


```R
head(tbl.all[,c(4,8,12,16,18)],10)
cor(tbl.all[c(4,8,12,16,18)])
```


<table>
<thead><tr><th></th><th scope=col>spearman.asinh</th><th scope=col>ridge.asinh</th><th scope=col>lpv.asinh</th><th scope=col>ensemble.asinh</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td> 0.9415501 </td><td> 0.04472347</td><td>48.38330171</td><td> 6.854995  </td><td> 0.9232394 </td></tr>
	<tr><th scope=row>2</th><td> 0.9205737 </td><td> 0.03628105</td><td>54.05992920</td><td>22.979982  </td><td> 0.9047599 </td></tr>
	<tr><th scope=row>3</th><td> 0.8903669 </td><td> 0.04830726</td><td> 0.01903403</td><td> 2.694855  </td><td> 0.8628270 </td></tr>
	<tr><th scope=row>4</th><td> 0.8765615 </td><td> 0.05485779</td><td>48.89733145</td><td> 2.407834  </td><td> 0.8423299 </td></tr>
	<tr><th scope=row>5</th><td> 0.8531967 </td><td> 0.07325669</td><td> 0.29125507</td><td> 2.019311  </td><td> 0.8278686 </td></tr>
	<tr><th scope=row>6</th><td> 0.8314088 </td><td> 0.05072036</td><td>29.16804671</td><td>16.120648  </td><td> 0.8196172 </td></tr>
	<tr><th scope=row>8</th><td> 0.8147850 </td><td> 0.04597392</td><td>37.06080200</td><td> 0.000000  </td><td> 0.7704362 </td></tr>
	<tr><th scope=row>7</th><td> 0.8169683 </td><td> 0.02333311</td><td> 0.06635882</td><td> 0.000000  </td><td> 0.7309450 </td></tr>
	<tr><th scope=row>10</th><td>-0.7672384 </td><td>-0.03368853</td><td>13.70820824</td><td> 0.000000  </td><td>-0.7273121 </td></tr>
	<tr><th scope=row>14</th><td> 0.7019364 </td><td> 0.04498817</td><td> 0.52651363</td><td>86.919009  </td><td> 0.7241906 </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>spearman.asinh</th><th scope=col>ridge.asinh</th><th scope=col>lpv.asinh</th><th scope=col>ensemble.asinh</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>spearman.asinh</th><td>1.0000000 </td><td>0.81964553</td><td>0.35552775</td><td>0.16380537</td><td>0.9949736 </td></tr>
	<tr><th scope=row>ridge.asinh</th><td>0.8196455 </td><td>1.00000000</td><td>0.24151081</td><td>0.03238692</td><td>0.8259104 </td></tr>
	<tr><th scope=row>lpv.asinh</th><td>0.3555278 </td><td>0.24151081</td><td>1.00000000</td><td>0.03731357</td><td>0.3807169 </td></tr>
	<tr><th scope=row>ensemble.asinh</th><td>0.1638054 </td><td>0.03238692</td><td>0.03731357</td><td>1.00000000</td><td>0.1809052 </td></tr>
	<tr><th scope=row>gene.cor</th><td>0.9949736 </td><td>0.82591042</td><td>0.38071693</td><td>0.18090515</td><td>1.0000000 </td></tr>
</tbody>
</table>



These results show about the ranking as the log2 results. Ridge, lassopv, and ensemble all correlated slightly less, while spearman continued to do so well. 

### voom Dataset


```R
head(tbl.all[,c(5,9,13,17,18)],10)
cor(tbl.all[c(5,9,13,17,18)])
```


<table>
<thead><tr><th></th><th scope=col>spearman.voom</th><th scope=col>ridge.voom</th><th scope=col>lpv.voom</th><th scope=col>ensemble.voom</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td> 0.9416976  </td><td> 0.04533243 </td><td>48.973631221</td><td> 6.967867   </td><td> 0.9232394  </td></tr>
	<tr><th scope=row>2</th><td> 0.9225219  </td><td> 0.03813049 </td><td>54.624312937</td><td>23.164242   </td><td> 0.9047599  </td></tr>
	<tr><th scope=row>3</th><td> 0.8931333  </td><td> 0.04796167 </td><td> 0.002191662</td><td> 2.469122   </td><td> 0.8628270  </td></tr>
	<tr><th scope=row>4</th><td> 0.8756181  </td><td> 0.05723240 </td><td>49.424139376</td><td> 2.393848   </td><td> 0.8423299  </td></tr>
	<tr><th scope=row>5</th><td> 0.8721602  </td><td> 0.07547096 </td><td> 0.135304746</td><td> 2.054066   </td><td> 0.8278686  </td></tr>
	<tr><th scope=row>6</th><td> 0.8365842  </td><td> 0.05313104 </td><td>29.826784860</td><td>16.252081   </td><td> 0.8196172  </td></tr>
	<tr><th scope=row>8</th><td> 0.8193281  </td><td> 0.04443798 </td><td>37.505220279</td><td> 0.000000   </td><td> 0.7704362  </td></tr>
	<tr><th scope=row>7</th><td> 0.8194466  </td><td> 0.02897746 </td><td> 0.016631006</td><td> 0.000000   </td><td> 0.7309450  </td></tr>
	<tr><th scope=row>10</th><td>-0.7615027  </td><td>-0.03612586 </td><td>15.936759724</td><td> 0.000000   </td><td>-0.7273121  </td></tr>
	<tr><th scope=row>14</th><td> 0.7219303  </td><td> 0.05327404 </td><td> 0.591778019</td><td>87.605358   </td><td> 0.7241906  </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>spearman.voom</th><th scope=col>ridge.voom</th><th scope=col>lpv.voom</th><th scope=col>ensemble.voom</th><th scope=col>gene.cor</th></tr></thead>
<tbody>
	<tr><th scope=row>spearman.voom</th><td>1.0000000 </td><td>0.8150102 </td><td>0.27120259</td><td>0.16672201</td><td>0.9856114 </td></tr>
	<tr><th scope=row>ridge.voom</th><td>0.8150102 </td><td>1.0000000 </td><td>0.16526521</td><td>0.10596919</td><td>0.8297440 </td></tr>
	<tr><th scope=row>lpv.voom</th><td>0.2712026 </td><td>0.1652652 </td><td>1.00000000</td><td>0.03108823</td><td>0.3492400 </td></tr>
	<tr><th scope=row>ensemble.voom</th><td>0.1667220 </td><td>0.1059692 </td><td>0.03108823</td><td>1.00000000</td><td>0.1843120 </td></tr>
	<tr><th scope=row>gene.cor</th><td>0.9856114 </td><td>0.8297440 </td><td>0.34923995</td><td>0.18431202</td><td>1.0000000 </td></tr>
</tbody>
</table>



Each solver correlated slightly better in the voom dataset, but the order of correlation has not changed. 

## Conclusions

After comparing the correlation of each solver with each dataset, and then each dataset with each solver, we have detailled the effectiveness of each solver (in comparison to the rest).

In terms of correlation, spearman clealy rose above the other solvers. Within each dataset, the spearman data remained similar and correlated to the given genetic correlation values closely. Ridge also performed noticeably well, with slightly less correlated data in both sets of tests. Ensemble's data correlated extremely closely when compared between datasets, but correlated the least out of all 4 solvers when compared to them. Lassopv did not correlate very closely in either set of tests.

These conclusions may be useful to anyone using TReNA and interested in the pros and cons of different solver usage. When making a decision between these solvers, take note of the fact that spearman and ensemble are more consistent when analyzing data across multiple transformations. Spearman and ridge seem to be the most closely correlated with gene correlations. 

These tests do NOT mean that spearman is the "best" solver, as this was a short examination of a small dataset in one particular case. However, it does provide a constant environment in which the different solvers are comparable.


```R

```
