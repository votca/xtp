/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSe-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/eigen.h>
#include <votca/xtp/gauss_hermite_quadrature_constants.h>

using namespace std;

namespace votca {
namespace xtp {

const Eigen::VectorXd& Gauss_Hermite_Quadrature_Constants::getPoints(
    Index order) {
  if (!this->_filled_Points) {
    this->FillPoints();
    _filled_Points = true;
  }
  if (_map_points.count(order) == 0){
    throw invalid_argument("Order " + std::to_string(order) +
                           " not in range {8,10,12,...,20}.");}
  return _map_points.at(order);
}

const Eigen::VectorXd& Gauss_Hermite_Quadrature_Constants::getAdaptedWeights(
    Index order) {
  if (!this->_filled_AdaptedWeights) {
    this->FillAdaptedWeights();
    _filled_AdaptedWeights = true;
  }
  if (_map_AdaptedWeights.count(order) == 0){
    throw invalid_argument("Order " + std::to_string(order) +
                           " not in range {8,10,12,...,20}.");}
  return _map_AdaptedWeights.at(order);
}

void Gauss_Hermite_Quadrature_Constants::FillPoints() {
  Eigen::VectorXd points_8(8);
  points_8 << -2.9306374202572440192235027052435991461994485855216,
      -1.9816567566958429258546306397693095686949116340382,
      -1.1571937124467801947207657790631002434520047595976,
      -0.38118699020732211685471888558369141763186003151966,
      0.38118699020732211685471888558369141763186003151966,
      1.1571937124467801947207657790631002434520047595976,
      1.9816567566958429258546306397693095686949116340382,
      2.9306374202572440192235027052435991461994485855216;
  _map_points[8] = points_8;
  Eigen::VectorXd points_10(10);
  points_10 << -3.4361591188377376033267254943191213848406783093902,
      -2.5327316742327897964089607977547934803078465081567,
      -1.7566836492998817734514012201061567632954744937389,
      -1.0366108297895136541774919167592090162982561106572,
      -0.34290132722370460878916502555725803120830265867773,
      0.34290132722370460878916502555725803120830265867773,
      1.0366108297895136541774919167592090162982561106572,
      1.7566836492998817734514012201061567632954744937389,
      2.5327316742327897964089607977547934803078465081567,
      3.4361591188377376033267254943191213848406783093902;
  _map_points[10] = points_10;
  Eigen::VectorXd points_12(12);
  points_12 << -3.8897248978697819192716427472441917601481810141665,
      -3.0206370251208897717106793751767636460038065983103,
      -2.279507080501059900187728569424341135039156593214,
      -1.5976826351526047967096627709045767082402869887687,
      -0.94778839124016374370457813106013651367915379591455,
      -0.31424037625435911127661163409533712897211083832271,
      0.3142403762543591112766116340953371289721108383227,
      0.94778839124016374370457813106013651367915379591455,
      1.5976826351526047967096627709045767082402869887687,
      2.279507080501059900187728569424341135039156593214,
      3.0206370251208897717106793751767636460038065983103,
      3.8897248978697819192716427472441917601481810141665;
  _map_points[12] = points_12;
  Eigen::VectorXd points_14(14);
  points_14 << -4.3044485704736318126212981003689427253218103879463,
      -3.4626569336022705502089173611504324421390765376615,
      -2.7484707249854025686249985241464309299737954499322,
      -2.0951832585077168157349727263032374554797833544483,
      -1.4766827311411408705835065442050785826659377346378,
      -0.87871378732939941611467931186079607260338557501491,
      -0.29174551067256207844611307579938114089108844043221,
      0.29174551067256207844611307579938114089108844043221,
      0.87871378732939941611467931186079607260338557501491,
      1.4766827311411408705835065442050785826659377346378,
      2.0951832585077168157349727263032374554797833544483,
      2.7484707249854025686249985241464309299737954499322,
      3.4626569336022705502089173611504324421390765376615,
      4.3044485704736318126212981003689427253218103879463;
  _map_points[14] = points_14;
  Eigen::VectorXd points_16(16);
  points_16 << -4.6887389393058183646884986487456108906323917275435,
      -3.8694479048601226987194240980148123970065994702497,
      -3.176999161979956026813994559263696476791667372742,
      -2.5462021578474813621593287054458941245194651994058,
      -1.9517879909162539774346554149598874927683260586599,
      -1.3802585391988807963720896696945820384406817699701,
      -0.82295144914465589258245449673394264076350817254127,
      -0.2734810461381524521582804019650150339298178803279,
      0.2734810461381524521582804019650150339298178803279,
      0.82295144914465589258245449673394264076350817254127,
      1.3802585391988807963720896696945820384406817699701,
      1.9517879909162539774346554149598874927683260586599,
      2.5462021578474813621593287054458941245194651994058,
      3.1769991619799560268139945592636964767916673727419,
      3.8694479048601226987194240980148123970065994702497,
      4.6887389393058183646884986487456108906323917275435;
  _map_points[16] = points_16;
  Eigen::VectorXd points_18(18);
  points_18 << -5.0483640088744667683720375788536521210964511478335,
      -4.2481178735681264630234201609020812454613481259571,
      -3.5737690684862660795006759937718894579393109766802,
      -2.9613775055316068447786325490618382146613988430531,
      -2.3862990891666860002645930142399451367070333328477,
      -1.8355316042616288922538394440906011492669936367257,
      -1.3009208583896173656662655543926105802181346396612,
      -0.7766829192674116613166594622838522947287132144905,
      -0.25826775051909675925811609871057963300171513532774,
      0.2582677505190967592581160987105796330017151353277,
      0.7766829192674116613166594622838522947287132144905,
      1.3009208583896173656662655543926105802181346396612,
      1.8355316042616288922538394440906011492669936367257,
      2.3862990891666860002645930142399451367070333328477,
      2.9613775055316068447786325490618382146613988430531,
      3.5737690684862660795006759937718894579393109766802,
      4.248117873568126463023420160902081245461348125957,
      5.0483640088744667683720375788536521210964511478335;
  _map_points[18] = points_18;
  Eigen::VectorXd points_20(20);
  points_20 << -5.3874808900112328620169004106811207539962864490659,
      -4.6036824495507442730776752489783475851133984877619,
      -3.944764040115625210375628800524411807149768127888,
      -3.3478545673832163269149245229964636985104785902937,
      -2.788806058428130480525033756403185410670698887902,
      -2.2549740020892755230823333447345651280822653160264,
      -1.7385377121165862067808656621364064429514094196004,
      -1.2340762153953230078858183469594102295854459300694,
      -0.7374737285453943587056051442521042290772162039768,
      -0.24534070830090124990383653063361662396613385130349,
      0.24534070830090124990383653063361662396613385130349,
      0.73747372854539435870560514425210422907721620397679,
      1.2340762153953230078858183469594102295854459300694,
      1.7385377121165862067808656621364064429514094196004,
      2.2549740020892755230823333447345651280822653160264,
      2.788806058428130480525033756403185410670698887902,
      3.3478545673832163269149245229964636985104785902937,
      3.9447640401156252103756288005244118071497681278883,
      4.603682449550744273077675248978347585113398487762,
      5.3874808900112328620169004106811207539962864490659;
  _map_points[20] = points_20;
  Eigen::VectorXd points_30(30);
  points_30 << -6.863345293529891581061,
-6.138279220123934620395,
-5.533147151567495725118,
-4.988918968589943944486,
-4.48305535709251834189,
-4.003908603861228815228,
-3.544443873155349886925,
-3.099970529586441748689,
-2.667132124535617200571,
-2.243391467761504072473,
-1.826741143603688038836,
-1.415527800198188511941,
-1.008338271046723461805,
-0.6039210586255523077782,
-0.2011285765488714855458,
0.2011285765488714855458,
0.6039210586255523077782,
1.008338271046723461805,
1.415527800198188511941,
1.826741143603688038836,
2.243391467761504072473,
2.667132124535617200571,
3.099970529586441748689,
3.544443873155349886925,
4.003908603861228815228,
4.48305535709251834189,
4.988918968589943944486,
5.533147151567495725118,
6.138279220123934620395,
6.863345293529891581061;
  _map_points[30] = points_30;

Eigen::VectorXd points_35(35);
  points_35 << -7.504021146448936080097,
-6.79960941328413015723,
-6.212973747633716873302,
-5.6864689480904415917,
-5.198099346197752624161,
-4.736518477413210800452,
-4.294895814492763231508,
-3.8687007309691543826,
-3.454716495751990827259,
-3.050538420430446690874,
-2.654292781197172035651,
-2.264467501042568649446,
-1.879803988730917078884,
-1.499224488611730168455,
-1.121780990720302616783,
-0.7466176398798670171184,
-0.3729417170496168611453,
0,
0.3729417170496168611453,
0.7466176398798670171184,
1.121780990720302616783,
1.499224488611730168455,
1.879803988730917078884,
2.264467501042568649446,
2.654292781197172035651,
3.050538420430446690874,
3.454716495751990827259,
3.8687007309691543826,
4.294895814492763231508,
4.73651847741321080045,
5.198099346197752624161,
5.686468948090441591697,
6.212973747633716873302,
6.79960941328413015723,
7.504021146448936080097;
_map_points[35] = points_35;

Eigen::VectorXd points_50(50);
  points_50 << -9.182406958129317366347,
-8.52277103091780418914,
-7.97562236820563655424,
-7.486409429864194266821,
-7.034323509770610648808,
-6.60864797385535900613,
-6.20295251927467161631,
-5.812994675420406059157,
-5.435786087224948141616,
-5.069117584917235032452,
-4.71129366616904278739,
-4.360973160454578664322,
-4.017068172858134387881,
-3.678677062515269281719,
-3.345038313937891090222,
-3.01549776957452241886,
-2.689484702267745072548,
-2.366493904298663828905,
-2.046071968686409207851,
-1.72780654751589855853,
-1.411317754898300062015,
-1.096251128957681642351,
-0.7822717295546068858116,
-0.4690590566782360862441,
-0.1563025468894686754381,
0.1563025468894686754381,
0.4690590566782360862441,
0.7822717295546068858116,
1.096251128957681642351,
1.411317754898300062015,
1.72780654751589855853,
2.046071968686409207851,
2.366493904298663828905,
2.689484702267745072548,
3.01549776957452241886,
3.345038313937891090222,
3.678677062515269281719,
4.017068172858134387881,
4.360973160454578664322,
4.711293666169042787394,
5.069117584917235032452,
5.435786087224948141616,
5.812994675420406059157,
6.202952519274671616315,
6.608647973855359006128,
7.034323509770610648808,
7.486409429864194266821,
7.97562236820563655424,
8.522771030917804189138,
9.182406958129317366347;
_map_points[50] = points_50;

  Eigen::VectorXd points_100(100);
  points_100 << -13.4064873381449101385
-12.82379974948780890634,
-12.3429642228596742951,
-11.91506194311416580199,
-11.52141540078703024169,
-11.1524043855851252649,
-10.80226075368471459482,
-10.46718542134281214178,
-10.14450994129284546989,
-9.83226980777796909436,
-9.5289658233901148047,
-9.233420890219161550478,
-8.9446892173254744788,
-8.661996168134517714376,
-8.38469694041626507509,
-8.112247311162791917212,
-7.844182384460821168792,
-7.580100807857488884286,
-7.31965282230453531633,
-7.062531060248865437466,
-6.80846335285879641448,
-6.55720703192153931598,
-6.308544361112135121639,
-6.06227883261430263867,
-5.81823213520351704736,
-5.576241649329924103303,
-5.33615836013836049728,
-5.097845105089136247,
-4.861175091791210210046,
-4.626030635787155773074,
-4.392302078682684016745,
-4.159886855131030540068,
-3.92868868342767097201,
-3.698616859318491939797,
-3.469585636418589169768,
-3.241513679631012950359,
-3.014323580331155516715,
-2.78794142398198931319,
-2.562296402372608025056,
-2.337320463906878505013,
-2.11294799637118795203,
-1.889115537427008371494,
-1.665761508741509469867,
-1.442825970215932787703,
-1.220250391218953058821,
-0.9979774360981052439241,
-0.7759507615401457819751,
-0.554114823591616988233,
-0.3324146923422318070459,
-0.110795872422439482888,
0.1107958724224394828876,
0.332414692342231807046,
0.554114823591616988233,
0.775950761540145781975,
0.9979774360981052439241,
1.220250391218953058821,
1.442825970215932787703,
1.665761508741509469867,
1.889115537427008371494,
2.11294799637118795203,
2.337320463906878505013,
2.562296402372608025056,
2.78794142398198931319,
3.014323580331155516715,
3.241513679631012950359,
3.469585636418589169768,
3.6986168593184919398,
3.92868868342767097201,
4.159886855131030540068,
4.392302078682684016745,
4.626030635787155773074,
4.861175091791210210046,
5.097845105089136247,
5.33615836013836049728,
5.576241649329924103303,
5.818232135203517047362,
6.06227883261430263867,
6.308544361112135121639,
6.55720703192153931598,
6.80846335285879641448,
7.062531060248865437466,
7.31965282230453531633,
7.580100807857488884286,
7.84418238446082116879,
8.112247311162791917212,
8.38469694041626507509,
8.661996168134517714376,
8.9446892173254744788,
9.233420890219161550478,
9.5289658233901148047,
9.83226980777796909436,
10.14450994129284546989,
10.46718542134281214178,
10.80226075368471459482,
11.1524043855851252649,
11.52141540078703024169,
11.91506194311416580199,
12.3429642228596742951,
12.82379974948780890634,
13.4064873381449101385;
  _map_points[100] = points_100;

}

void Gauss_Hermite_Quadrature_Constants::FillAdaptedWeights() {
  Eigen::VectorXd AdaptedWeights_8(8);
  AdaptedWeights_8 << 1.0719301442479797564640154828686192710583035011087,
      0.8667526065633812222177945178444721007823603317873,
      0.7928900483864012509056467065099944670699652191022,
      0.7645441286517291990713914850625153797906970121608,
      0.7645441286517291990713914850625153797906970121608,
      0.7928900483864012509056467065099944670699652191022,
      0.8667526065633812222177945178444721007823603317873,
      1.0719301442479797564640154828686192710583035011087;
  _map_AdaptedWeights[8] = AdaptedWeights_8;
  Eigen::VectorXd AdaptedWeights_10(10);
  AdaptedWeights_10 << 1.0254516913657372330206021439998001007206764041771,
      0.820666126404816614571910879903882294986432537887,
      0.7414419319435649700820477018891891603635366640119,
      0.70329632310490617009834172682079856333223356450024,
      0.68708185395127336268656718422692784453580030087846,
      0.6870818539512733626865671842269278445358003008785,
      0.7032963231049061700983417268207985633322335645002,
      0.7414419319435649700820477018891891603635366640119,
      0.8206661264048166145719108799038822949864325378872,
      1.0254516913657372330206021439998001007206764041771;
  _map_AdaptedWeights[10] = AdaptedWeights_10;
  Eigen::VectorXd AdaptedWeights_12(12);
  AdaptedWeights_12 << 0.989699047092298099365013260061309904979902995276,
      0.7866439394633224644915380743541227886372855230251,
      0.7052203661122197557556715481998492923710656431374,
      0.662662773266871319239255150461100359348980754644,
      0.63962123202025660066934991623180131877302719880201,
      0.62930787436949282103571883185136329069537833837549,
      0.62930787436949282103571883185136329069537833837549,
      0.639621232020256600669349916231801318773027198802,
      0.66266277326687131923925515046110035934898075464401,
      0.70522036611221975575567154819984929237106564313742,
      0.7866439394633224644915380743541227886372855230251,
      0.9896990470922980993650132600613099049799029952755;
  _map_AdaptedWeights[12] = AdaptedWeights_12;
  Eigen::VectorXd AdaptedWeights_14(14);
  AdaptedWeights_14 << 0.960878703025659256530957856273245078574670135086,
      0.7599870873975660800632761513147237493826116686251,
      0.677706759192397010397799872260068292400349536379,
      0.63290060647233309082923787100002824017391731693,
      0.6063797391260974872052603343591461156764643412857,
      0.59110666704316242483025610616840494300223021601381,
      0.58406169052199623900075832593066558843366126221565,
      0.5840616905219962390007583259306655884336612622156,
      0.59110666704316242483025610616840494300223021601381,
      0.6063797391260974872052603343591461156764643412857,
      0.6329006064723330908292378710000282401739173169296,
      0.677706759192397010397799872260068292400349536379,
      0.75998708739756608006327615131472374938261166862506,
      0.960878703025659256530957856273245078574670135086;
  _map_AdaptedWeights[14] = AdaptedWeights_14;
  Eigen::VectorXd AdaptedWeights_16(16);
  AdaptedWeights_16 << 0.936874492884069357470004849342561171325068702881,
      0.738245622277681359882790450735840033408093284309,
      0.6557556728761177064705989221774875270966151707271,
      0.6097369582559972856242997148163131018517144433067,
      0.58124727540086389202624422184292731510412746076873,
      0.5632178290881998377123686952648505613118715089659,
      0.55244195736745939041573957443504654808920309409138,
      0.54737520503784399928196163481039041583828307915463,
      0.54737520503784399928196163481039041583828307915463,
      0.55244195736745939041573957443504654808920309409138,
      0.5632178290881998377123686952648505613118715089659,
      0.5812472754008638920262442218429273151041274607687,
      0.6097369582559972856242997148163131018517144433067,
      0.65575567287611770647059892217748752709661517072711,
      0.738245622277681359882790450735840033408093284309,
      0.936874492884069357470004849342561171325068702881;
  _map_AdaptedWeights[16] = AdaptedWeights_16;
  Eigen::VectorXd AdaptedWeights_18(18);
  AdaptedWeights_18 << 0.9163935375519155947662880291953182438433621949995,
      0.719993383105314144400021542382132236413335255132,
      0.637630172006160187495403399833138046386032910798,
      0.5909530034631081328753352419791538626101162669642,
      0.561279045549804538746093740280135250066751747336,
      0.54157867866213244558778798431831136530920994802962,
      0.5285894429188008800743612444401935861230623382277,
      0.5206349466760646162215527041692833691475359172075,
      0.5168458364816212822612399422473770735127756668207,
      0.5168458364816212822612399422473770735127756668207,
      0.52063494667606461622155270416928336914753591720746,
      0.5285894429188008800743612444401935861230623382277,
      0.5415786786621324455877879843183113653092099480296,
      0.561279045549804538746093740280135250066751747336,
      0.590953003463108132875335241979153862610116266964,
      0.637630172006160187495403399833138046386032910798,
      0.7199933831053141444000215423821322364133352551319,
      0.916393537551915594766288029195318243843362194999;
  _map_AdaptedWeights[18] = AdaptedWeights_18;
  Eigen::VectorXd AdaptedWeights_20(20);
  AdaptedWeights_20 << 0.898591961453191416420545454828522781943430553196,
      0.7043329611769424087138204002833519583892017111602,
      0.622278696191412280649317396569296403505856898279,
      0.575262442852503182139156253731619647430862573735,
      0.544851742364520005199359051276235186443829263999,
      0.524080350948557612676557263766842812596750345434,
      0.50967902711745800667286665991982639104787118706115,
      0.4999208713362905172649779639324870322555667062887,
      0.49384338527205292781460899199861133161984929633202,
      0.49092150066674582427959429560835390391359129040998,
      0.49092150066674582427959429560835390391359129040998,
      0.493843385272052927814608991998611331619849296332,
      0.4999208713362905172649779639324870322555667062887,
      0.5096790271174580066728666599198263910478711870611,
      0.5240803509485576126765572637668428125967503454339,
      0.544851742364520005199359051276235186443829263999,
      0.575262442852503182139156253731619647430862573735,
      0.622278696191412280649317396569296403505856898279,
      0.70433296117694240871382040028335195838920171116,
      0.898591961453191416420545454828522781943430553196;
  _map_AdaptedWeights[20] = AdaptedWeights_20;
  
 Eigen::VectorXd AdaptedWeights_30(30);
  AdaptedWeights_30 << 0.83424747101276179534,
0.64909798155426670071,
0.56940269194964050397,
0.52252568933135454964,
0.491057995832882696506,
0.46837481256472881677,
0.45132103599118862129,
0.438177022652683703695,
0.4279180629327437485828,
0.4198950037368240886418,
0.413679363611138937184,
0.4089815750035316024972,
0.4056051233256844363121,
0.403419816924804022553,
0.402346066701902927115,
0.4023460667019029271154,
0.4034198169248040225528,
0.4056051233256844363121,
0.4089815750035316024972,
0.413679363611138937184,
0.4198950037368240886418,
0.427918062932743748583,
0.4381770226526837037,
0.45132103599118862129,
0.46837481256472881677,
0.4910579958328826965056,
0.52252568933135454964,
0.56940269194964050397,
0.64909798155426670071,
0.83424747101276179534;
_map_AdaptedWeights[30] = AdaptedWeights_30;

 Eigen::VectorXd AdaptedWeights_35(35);
  AdaptedWeights_35 << 0.81136071084102983378,
0.62989691645498025901,
0.551417853029805920067,
0.50497977495490646718,
0.47356764273572666754,
0.45070011483121789192,
0.433285296641151855372,
0.419633160917753598828,
0.40873163584413443683,
0.39993408957772014563,
0.39280745555744834385,
0.387051690293202808822,
0.382454126977783308389,
0.37886225962583710092,
0.3761669023632347819562,
0.3742915302886360867873,
0.373185514995801580782,
0.3728199731907248355968,
0.3731855149958015807822,
0.3742915302886360867873,
0.376166902363234781956,
0.37886225962583710092,
0.382454126977783308389,
0.387051690293202808822,
0.39280745555744834385,
0.39993408957772014563,
0.408731635844134436828,
0.419633160917753598828,
0.43328529664115185537,
0.4507001148312178919154,
0.473567642735726667542,
0.50497977495490646718,
0.551417853029805920067,
0.62989691645498025901,
0.81136071084102983378;
_map_AdaptedWeights[35] = AdaptedWeights_35;

 Eigen::VectorXd AdaptedWeights_50(50);
  AdaptedWeights_50 << 0.7613486911180767503,
0.5886052973772898403,
0.51330479785154088953,
0.46832621194254749301,
0.43755328230004755951,
0.41483882105902249664,
0.39724494377584631487,
0.38316139219647935478,
0.37161977124991502943,
0.36199726509161745079,
0.35387246997946172103,
0.34694876959937838205,
0.341010727258747850608,
0.335897876837685552073,
0.331488254247521860023,
0.327687661354606419534,
0.3244224487088668268899,
0.32163453796205597951,
0.319277915993180742134,
0.317316124349604325531,
0.315720440230093698585,
0.3144685508824719215511,
0.3135435899896416108314,
0.312933448052389891275,
0.3126302980303591226471,
0.3126302980303591226471,
0.312933448052389891275,
0.313543589989641610831,
0.314468550882471921551,
0.315720440230093698585,
0.317316124349604325531,
0.319277915993180742134,
0.321634537962055979507,
0.32442244870886682689,
0.32768766135460641953,
0.331488254247521860023,
0.335897876837685552073,
0.34101072725874785061,
0.346948769599378382046,
0.353872469979461721026,
0.36199726509161745079,
0.371619771249915029427,
0.38316139219647935478,
0.39724494377584631487,
0.41483882105902249664,
0.43755328230004755951,
0.46832621194254749301,
0.51330479785154088953,
0.5886052973772898403,
0.7613486911180767503;
_map_AdaptedWeights[50] = AdaptedWeights_50;

  Eigen::VectorXd AdaptedWeights_100(100);
  AdaptedWeights_100 << 0.674353552420908144,
0.5185068072702001877,
0.449993171054416591,
0.4086886584419015471,
0.380137438800575956,
0.3588184090711120834,
0.3420892623142471455,
0.3284994846454574786,
0.31717511097185385938,
0.30755272850167372943,
0.29924995803907909358,
0.29199656396027608226,
0.28559521453136788494,
0.27989787254508465858,
0.274790938326829146069,
0.2701855435339907877,
0.266011005285956696413,
0.26221028944375895952,
0.25873679090438086577,
0.25555200056201943528,
0.25262378339141757328,
0.24992508660128928222,
0.24743295612683529342,
0.245127777913567834895,
0.24299268557910773048,
0.24101309292238722418,
0.23917632129955903596,
0.2374712999204616475105,
0.23588832279469925466,
0.23441885012175822832,
0.233055344869710310074,
0.2317911374537038421,
0.230620313034512848758,
0.22953761716486989448,
0.228538376426246803542,
0.22761843139840206564,
0.226774079843685723664,
0.22600202840780008501,
0.225299351467708196188,
0.224663456017245197001,
0.224092051687791109734,
0.223583125167177609196,
0.2231349184139684965647,
0.222745910173357235135,
0.2224148003905443508998,
0.222140497191773490166,
0.221922106165498110136,
0.221758921729023785805,
0.2216504204115155471103,
0.2215962559241832702211,
0.221596255924183270221,
0.22165042041151554711,
0.2217589217290237858046,
0.221922106165498110136,
0.2221404971917734901661,
0.2224148003905443508998,
0.222745910173357235135,
0.223134918413968496565,
0.223583125167177609196,
0.224092051687791109734,
0.2246634560172451970005,
0.225299351467708196188,
0.226002028407800085009,
0.226774079843685723664,
0.22761843139840206564,
0.22853837642624680354,
0.22953761716486989448,
0.23062031303451284876,
0.231791137453703842101,
0.233055344869710310074,
0.23441885012175822832,
0.23588832279469925466,
0.23747129992046164751,
0.23917632129955903596,
0.24101309292238722418,
0.24299268557910773048,
0.245127777913567834895,
0.24743295612683529342,
0.24992508660128928222,
0.25262378339141757328,
0.25555200056201943528,
0.25873679090438086577,
0.26221028944375895952,
0.26601100528595669641,
0.270185543533990787703,
0.27479093832682914607,
0.27989787254508465858,
0.28559521453136788494,
0.29199656396027608226,
0.29924995803907909358,
0.30755272850167372943,
0.3171751109718538594,
0.32849948464545747856,
0.34208926231424714546,
0.3588184090711120834,
0.380137438800575956,
0.4086886584419015471,
0.4499931710544165908,
0.51850680727020018773,
0.674353552420908144;
  _map_AdaptedWeights[100] = AdaptedWeights_100;
}
}  // namespace xtp
}  // namespace votca
