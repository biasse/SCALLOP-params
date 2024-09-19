#pragma once

#include <NTL/ZZ.h>
#include <vector>

#include "fp2.hpp"

const NTL::ZZ p = NTL::conv<NTL::ZZ>("23771452727449498898004557547225354988120021364701398066313318892838230235651289792849929330578203612228877729959916570183072595360843544888617115522542149031138111127294471888062615222219308670025415410876727090874823170475481539518868548152615795234974491096957906735993259268864790983747139027153544490397993988855511337575848418797409297084369809840481360670399419714736596634148006572830334348569025487101216814643806207");
const NTL::ZZ lampos = NTL::conv<NTL::ZZ>("1205794941487557338817933421586945351024768814661670785658311427689589012033869192016273739935809183967388006252240405091621900299263526317091679942798365976676737629375085072508815266329627853");
const NTL::ZZ lamneg = NTL::conv<NTL::ZZ>("1412970443851765635529985925689555099829072726553431131483688314861114903531396612132134034171272046947455092273733495645106478991766762985675934573657981153961144064626711799861572492875849821");

const NTL::ZZ A_re = NTL::conv<NTL::ZZ>("20775255613076975120462396342329376411570862479699901766707191690065658286185162486429557613180766319688053566597747061294388377457672235013594851637510152387255470818124717477325816223876908965954280785355367308691830276201776524763446492381170644054738990196592456723505609278220709190415549424423999745780553298197610684071007828244690481586336035413382617664304485043691847158933172950081774757002971874877902784836075176");
const NTL::ZZ A_im = NTL::conv<NTL::ZZ>("19541278133547132606637252739398665909927598778915435471555577348188176951288345329534809443256304824951352423878253112121747483455127493560936729568065499127454131191819031846814596228456554039607836917905720595469634405166456665799482879701144702089243016512225242720226837904201189802140210141269755522455665370837872676048031947362754539193227137077101983952989554925309061319448701760525965937491041677194889871288099236");

const NTL::ZZ P_re =  NTL::conv<NTL::ZZ>("23371311189665093790301003936368030246153473518710833041304373495431606658495101574593217419622474036025807417717039900273368534584417232518839849151657204941782006390815098581728512188879873404780806252911043417496367127541522243078523849652624382740321805812750608576678185095431933412252897204214683745234372744610165562198002175517153964231175093080740432839358774915606007657620225929404601286352946231444781059911742519");
const NTL::ZZ P_im =  NTL::conv<NTL::ZZ>("7563428169371363999191281502836763896000439770840311627947928503066354991902243896473979982851197429599803683292013163047024540955370030978812343403695384745302324719650130692143731555836759647018568752084517619245790079761297337246064274844868114132783079055750703645597032582472229140541281131358628004889147884545983848796538056057321579429277340106093330431849419894451663767056370554191519369122239394491669897494106104");
const NTL::ZZ Q_re =  NTL::conv<NTL::ZZ>("18871273384884259058789600913747058394000233620320688953523699028502191732875545477006532663148348586031792362513472957611543537543288865259429693273083488817321540059963998588006856495178488022934394166711785991415452267798693791536352261214364424640007305665120479903777074411937038249352704124012559023338957143317491921632832369561725784976444737953041723405483953902730479789132935635322252535550131008981835762694045494");
const NTL::ZZ Q_im =  NTL::conv<NTL::ZZ>("12718636333070363731375360665002669496491031334512381380936867996170191305422150183121018637263856678122032888327186929335065029159931118999996184704312527104089161670494122932252375731265787556806632506582792415931861793839201548735127625807184558744519890938082594978364813403131059055417893742387458288013060604293627826061699363963215707166775163684803327031220058916344375989477818935471972490025876724159639192791926515");
const NTL::ZZ Qm_re =  NTL::conv<NTL::ZZ>("19092161908971671246838211110203725193208266607734394973211344168464751049823293032791588722427684719251132136739429047382397621995428874160693071923230642192423345034797870568940777958848536897335582452043260500923035330605035895897283932979583852296568385474481127097406305301846587286735935595270809460969987312404283139077924070952662408836757296799493675696816999892539403085265338608649099682375992308143657951270004466");
const NTL::ZZ Qm_im =  NTL::conv<NTL::ZZ>("16855050570179969570900534716296447021332044718563453081397527900729372734483817168806948227102649750652281952173634425264329568195332701717427639219274377210939897300360396062312274904687250992365011882198055867133282424489587918023264936025269560777272275967457832286072778984881023215675066295412458055878413212388255626427292288526686595570367655017644365178177030984353842509422246234431917017601688111753506780015326222");

const std::vector<int> ells{11, 23, 37, 47, 53, 59, 67, 73, 97, 103, 127, 131, 137, 139, 179, 193, 197, 223, 229, 233, 257, 269, 281, 283, 293, 313, 331, 347, 349, 353, 373, 383, 397, 419, 431, 443, 479, 491, 563, 577, 593, 601, 607, 613, 641, 647, 653, 659, 691, 727, 743, 751, 757, 761, 811, 821, 823, 829, 839, 857, 859, 877, 883, 887, 953, 983, 991, 997, 1031, 1061, 1063, 1069, 1087, 1093, 1097};
// 2 isog -> const std::vector<int> strat{256, 128, 64, 32, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 128, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1};
const std::vector<int> strat{157, 103, 58, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 45, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 18, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 4, 2, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 71, 38, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1};
const int dim = 75;
const int N_ball = 20;
const int e = 762;