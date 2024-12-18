proof.all(False)

p = 4114534055211428858209763614504410563117701662267362937138137095586936058053761704842156133975358905714291000507768794648473124720678019767543502087984105654039860176528637932566910494321588571084727940922746613867782136100904100899254907911401022219784668284602450834434561264189531655905496525890573382474959042286886123544723208330612164409523725883983605918613518273414945403217222161064894795280154916448207951095491621077495805609660640803620150293725576124452440320886343175508418100084863559079428799766681877685071383771673711806607518344260958541082058104638843060207642392140545208993133008939248313207931109695817894592511
e = 762
B = QuaternionAlgebra(-1, -p)
i,j,k = B.gens()
O0 = B.maximal_order()

omega = 1832827684332475215473523438232831309177774886978512535974698570911369604699151107875446121559074470856162325880745766383242293239553824886289089398891540880222569935499346998589048849203406574145521169540424917988523539271567833574172297078821253610246946076670772916648107040797719390100519790897627262641282454683649730371689854402668571027616943858763342329950228548282868159012438018581416850038531760014452514324869538461260930969465953785455634366396271982746348423500078414182266292377619470577014647990002429998847200610475363634076413706320343626150729214765910589216653630468185680660416651415730334493810319468061039986231/3703584420421317772483294103983984255788644734551181071445961422091673356197739035533621198296498218136732562553643079540487416747623473599096681055653087036614392281291272719588162948853126099350338685863004191436628396431348500141676038247223060499481736762111223449615927880771118388904338240885264262053314220103524721239065994213298001420173261209173635603910579337964619493127940829684644326179340588525*i + 577318789112374288731069338708962744800035755084063208244645707003475694872724830357211154684448651819240797511899715514232221175608208266997793081271105980251679704130151626737974805502813285894624623687640151408452377622431233690685912573035187136796161775531716097653740352167242352102278190884840433045574550224/2664244134028040843246810572124292774432935547084635262930092016004288415350573039006077372674548231277687661851672626308972918449891608055420684203196651075503103200898400424994739936460377581977746091389370883994274107267533002957093258466595061027303464238444395610996799441173318592672605951114025223995932845531692683540954756873606133499200609168340016706563292852970756562716352647447026974311575*j + 558646372121734059507028928074482530274313966907240927861306687817085216423327801840711268784017043659326050790375375624524784612660645296275727620261206673128649618230276107943299097161651841251933674882191827700428401688656183423711867926885794574409063334599723688079849269769992909917961291100827337384114586430583/6441016383341422213014424528667798705719382147045532298166889429724649315126502670493254257906953422846491413136770573113891159561084301911472488792440151368025030054419604729718544258875001911913632497153050767715875472054519130681175718690822713912142150890628214694984222401341075458964066505887416107918807339310477776067940859501387828556823062972475888006801007544286294770657288399451555349877114067*k
path = 3703584420421317772483294103983984255788644734551181071445961422091673356197739035533621198296498218136732562553643079540487416747623473599096681055653087036614392281291272719588162948853126099350338685863004191436628396431348500141676038247223060499481736762111223449615927880771118388904338240885264262053314220103524721239065994213298001420173261209173635603910579337964619493127940829684644326179340588525

assert omega*path in O0

I = O0*path + O0*(omega*path)

lamsqr = 2^(e*2) - omega.reduced_norm()

lam = lamsqr.sqrt()
w = omega - lam
assert w.reduced_norm() == 2**(e*2)

O = I.right_order()
assert w in O

print(f"w = {lam}")