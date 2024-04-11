#pragma once

#include <NTL/ZZ.h>
#include <vector>

#include "fp2.hpp"

const NTL::ZZ p = NTL::conv<NTL::ZZ>("167286323857221689112346016933207258999493176647479781908348180838625562682489086996433613891517156716513168813389511283523303877305870043625319556074092350622078682027926127121250390130545025689333161800187261717666012808114506310888652381411304126074312458239");
const NTL::ZZ lampos = NTL::conv<NTL::ZZ>("6371085725815258476752908878975343421013800487502345682680163749521241114919138170198188871492673549970978336280604153635105449368844094666212882510222463365349494215188111078002567");
const NTL::ZZ lamneg = NTL::conv<NTL::ZZ>("53825284682229889663454748309833518725356896733497612930262010008065229280626999849480050081613790224330260956400492502091056737710284118850381855605828771970300358158725245400615418");

const NTL::ZZ A_re = NTL::conv<NTL::ZZ>("96177889592231928148956246270625562523137886102010052670759059710202697191962387742141410331051700000446817311207205914182213516215162823927337435936179195242608216981927342398422767478475483417343265983445206684917540848268112973254905224910041523422577643332");
const NTL::ZZ A_im = NTL::conv<NTL::ZZ>("23685146305789021641363358835934420664343939425256747386805118608062199963728183817569500571913148698343968280079924384030333496851055862339208650812640352534420861071335964235423464820610991477429182641790598375349536584897532814977077649279996871378365350617");

const NTL::ZZ P_re =  NTL::conv<NTL::ZZ>("41822696800159218449979642297569007521287323847135525637805187117401212315840417805686273337005043641885142255827804970573398737826568047696956788719156374222467923715724461827187197774216359457772241568067451410104393376266142755667327405827610575709481444282");
const NTL::ZZ P_im =  NTL::conv<NTL::ZZ>("140022157426865414875023043049491984408056107891697367465507060466455830737816622565980558919564319223298887797470322871872472813207458680733461393608737835761074827996009541217073565088662561428011367036365767084667257519780084226184466300948104098553534241659");
const NTL::ZZ Q_re =  NTL::conv<NTL::ZZ>("45149074969443068493508830065578203608147116007752997886465094565740885315104236745340796000591630463564619657739144155041327463053244034002667524941683842707244077607068412454083668215878740632192772542799956820910336652782400126367671237129488778903237285134");
const NTL::ZZ Q_im =  NTL::conv<NTL::ZZ>("35290597483765161135433356693008876944304907007647127168152658854489941337683286877738514808400004350710891196246586780231889658667859124342440644690616214348644234742690658734193704372692355954694674374801557259251813261050509480518408698700004886392218485352");
const NTL::ZZ Qm_re =  NTL::conv<NTL::ZZ>("149836088380760014816061204171737844102620848911066437338891731226783314076359613985257797912965073507795900226577660578644310234731099424983501600349514308711897300673971439744727298478422694312539036058601712514898451621724651645723381600839241056681024730315");
const NTL::ZZ Qm_im =  NTL::conv<NTL::ZZ>("105646970810839919934582140110363867813700460914277269990920691102390081673011973687828273997512456646314443196942384551543752428277936048592842173616730276452053610737352086669858929956534385749973130309753076738896246677235094437139144875827335095126803328813");

const std::vector<int> ells{5, 7, 11, 17, 19, 23, 31, 41, 53, 61, 71, 79, 101, 109, 113, 131, 137, 149, 157, 163, 181, 193, 199, 223, 227, 229, 239, 251, 257, 269, 271, 283, 313, 331, 337, 349, 373, 419, 431, 443, 449, 457, 463, 479, 491, 503, 523, 557, 563, 571, 587, 593, 601, 613, 643, 647, 653, 683, 739, 743, 751, 769, 773, 797, 811, 839, 853, 857, 863, 877, 881, 883, 911, 937, 941};
const std::vector<int> strat{128, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1};

const int dim = 75;
const int N_ball = 7;