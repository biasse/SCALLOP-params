/* relations-1500.gp */

suffix = "1500";

default(parisize,"1G");

d = -4657663250639184493417860902365641674114776564696771899003406891564171427055
d0 = d \ 9;
factd = factor(d)[,1];

pol = x^2-d0;
nf = nfinit(pol);

addprimes([p | p <- factd, p>10^7]);
a=155853646279738557685783060518469225577;
N=2^126;
fsmall = 3^2 * 32
p1 = a;         \\ p1,p2,p3
p2 = a-N;       \\ have less than 128 bits
p3 = (a+N)\3;   \\ we can use pari
p4 = a^2-3*N^2; \\ p4 has 250 bits, use cado
flarge = p1*p2*p3*p4;
f = fsmall * flarge;
addprimes([p1,p2,p3,p4]);

m0small = znstar(fsmall).cyc[1];
m1small = 2^3 * 59;
m1a = 120847206871530317;
m1b = 2732362722177477299;
m2small = 2^3 * 3^2 * 43 * 234539;
m2 = 97479496874890883339561132173;
m3small = 2 * 3 * 11 * 43;
m3a = 1809194370466573;
m3b = 15640895719669203679;
m4small = 2^4 * 3^2 * 5 * 13 * 43;
m4a = 6225897700316611;
m4b = 1029349126277242586407648729407215118817378323925924073;
Msmall = lcm([m0small,m1small,m2small,m3small,m4small]);
Mlarge = [m1a,m1b,m2,m3a,m3b,m4a,m4b];
addprimes(Mlarge);

\\om = ***
cc = Mat([nfalgtobasis(nf,subst(b,'x,-'x)) | b <- nf.zk]); \\matrix of complex conjugation

\\primes to avoid in factorbasis = 2,3 ?????

/* checks */
\\if(nfeltnorm(nf,om)!=2^***, error("wrong omega norm!"));
if(nf.zk != [1,(x-1)/2], error("ZK basis has changed!"));
if(idealstar(nf,f,0).cyc[1] != Msmall * factorback(Mlarge), error("wrong decomposition of Phi(f*ZK)"));
x1 = vectorv(2,i,random(10^10));
x2 = vectorv(2,i,random(10^10));
if(cc*nfeltmul(nf,x1,x2) != nfeltmul(nf,cc*x1,cc*x2), error("wrong automorphism"));


printtime(t) = 
{
  if(t<1000,
    print("[", t, " ms]");
  ,\\else
    print("[", t, " ms = ", strtime(t) "]");
  );
};

writematrix(filename, M) = \\fplll format (check?)
{
  my([m,n] = matsize(M), F);
  F = fileopen(filename,"w");
  for(i=1,m,
    filewrite1(F,"[");
    for(j=1,n, filewrite1(F,strprintf(" %d", M[i,j])));
    filewrite(F," ]");
  );
  fileclose(F);
};

d0 = - randomprime(2^128, Mod(3,4)); \\only for testing
system(strprintf("touch STRUCTURE-%s",suffix));
data = readvec(strprintf("STRUCTURE-%s",suffix));
{if(#data,
  print("loading already computed relations for maximal order.");
  h = data[1];
  cyc = data[2];
,\\else
  print("computing relations for maximal order.");
  res = quadclassunit(d0);
  h = res.no;
  cyc = res.cyc;
  write(strprintf("STRUCTURE-%s",suffix),h);
  write(strprintf("STRUCTURE-%s",suffix),cyc);
  system(strprintf("mv BASIS BASIS-%s",suffix));
  system(strprintf("mv REL REL-%s",suffix));
)};
FB = read(strprintf("BASIS-%s",suffix));
FB = Vec(FB);
M = readvec(strprintf("REL-%s",suffix));
M = matconcat(M);
print("h=",h);
print("cyc=",cyc);
if(h != vecprod(cyc), error("wrong structure"));
if(#M~ != #FB, error("dimension of relation matrix != size of factor base"));

/*
From Bill on the 1500 bits example:
After dividing d by 3^2, the class number is 18956989197469414174124120710182322176
and the class group SNF is [4628171190788431194854521657759356,4,2,2,2,2,2,2,2,2,2,2]
with generators:
[Qfb(44587,10449,2901732225536083321133626955518699816659799047490656157252741140021602),Qfb(8726019445724791095377793482388102477,-8442782027599716148886722981310291285,16869052120149269483909349896657032140),Qfb(12940188378841856893479649465,12940188378841856893479649465,9998272896206228020391280198069869879723608842),Qfb(63329583299744487007762391455,63329583299744487007762391455,2042955724619450476824004665217599844183472706),Qfb(1368263700914586557297419808623508087,1368263700914586557297419808623508087,94899521957641381742870930650481127918),Qfb(35966514436255,35966514436255,3597221937346257431631919880506880963078076056610969123857346),Qfb(13306246495908188253535,13306246495908188253535,9723217947282347764225940108005617122496406653156258),Qfb(6710175606588562069,6710175606588562069,19281095209034865569890642455529287993018511782150831756),Qfb(871023985,871023985,148537281370015714365644047179380793659098216601309648302759010598),Qfb(42468881890391952731892467472114455,42468881890391952731892467472114455,3046465550361715245988422818093952900606),Qfb(39218881,39218881,3298909388566628074864885233867643208010204578005830561137844645294),Qfb(4636809081749,4636809081749,27902709052488205726724282822891534293702937333530416500715276)],1]
*/


system(strprintf("touch LLLREL-%s",suffix));
Mlll = read(strprintf("LLLREL-%s",suffix));
{if(type(Mlll)=="t_INT",
  print("LLL reducing relations for maximal order.");
  H = mathnfmodid(M,2*h);
  Mlll = qflll(H,3); \\in place, no flatter
  write(strprintf("LLLREL-%s",suffix), Mlll);
,\\else
  print("loading LLL reduced relations for maximal order.");
)};
\\here we could check that mathnfmodid(Mlll,2*h) is consistent.


nb = 75; \\number of primes we want to keep
forbidden = [2,3];
deleteprimes = List();
keepprimes = List();
for(i=1,#FB,
  my(p = FB[i]);
  if(vecsearch(forbidden,p),
    listput(~deleteprimes,i);
  ,/* else if */#keepprimes<nb,
    listput(~keepprimes,i);
  ,\\else
    listput(~deleteprimes,i);
  );
);
keepprimes = Vec(keepprimes);
deleteprimes = Vec(deleteprimes);
subFB = vecextract(FB,keepprimes);
system(strprintf("touch SUBREL-%s-%s",suffix,nb));
Msub = read(strprintf("SUBREL-%s-%s",suffix,nb));
{if(type(Msub)=="t_INT",
  print("computing maximal order relations for a smaller factor base (", nb, " primes).");
  Mdel = vecextract(Mlll,deleteprimes,[1..#Mlll]);
  K = matkerint(Mdel);
  Msub = vecextract(Mlll,keepprimes,[1..#Mlll]);
  Msub = Msub*K;
  Msub = qflll(Msub,3);
  write(strprintf("SUBREL-%s-%s",suffix,nb),Msub);
,\\else
  print("loading maximal order relations for a smaller factor base (", nb, " primes).");
)};
\\could also store/load subFB, or check compatibility.


/*
  convention:
  D = -d
  chosen qf of norm p:
  p*x^2 + b*x*y + c*y^2
  with
  b = +- sqrt(D mod p) in [1,p-1]
  s.t. b=D mod 2 (here odd)
  Not sure about even primes.
*/
dec2 = idealprimedec(nf,2);
orderdec2 = [dec2[2],dec2[1]];
orderdec(p) =
{
  my(dec,b);
  if(p==2,return(orderdec2));
  dec = idealprimedec(nf,p);
  b = nfelttrace(nf,dec[1].gen[2])/2;
  if((b-d0)%2,
    dec
  ,\\else
    [dec[2],dec[1]]
  );
};

G = [2,-1;-1,(1-d0)/2]; \\det == -d0
reconstruct(C) =
{
  my(id,Gid,T,res);
  id = idealfactorback(nf, vector(nb,i,orderdec(B[i])[if(sign(C[i])>0,1,2)]), abs(C));
  Gid = id~*G*id;
  T = qflllgram(Gid);
  res = id*T[,1];
  if(nfeltnorm(nf,res) != factorback(B,abs(C)),
    error("wrong norm in reconstruct: ", res)
  );
  res
};

/*
L = Vec(M);
system(strprintf("mv ALGREL-ex2-%d ALGREL-ex2-%d-prev",nb,nb));
print1("reconstruct:");
t = getabstime();
{for(i=1,nb,
  if(i%20==0,print1(" ",round(100*i/nb),"%"));
  C = L[i];
  g = reconstruct(C);
  write(strprintf("ALGREL-ex2-%d",nb), g);
)};
t = getabstime()-t;
print(" done.");
printtime(t);
\\*/

/*
R = readvec(strprintf("ALGREL-ex2-%d",nb));
RQ = [[r,1;cc*r,-1] | r <- R];
*/

\\compute discrete logs at small primes with idealstar with cycmod
/*
bidsmall = idealstar(nf,f,1,Msmall);
bidcyc = [gcd(d,Msmall) | d <- bidsmall.cyc]; \\should be done by idealstar (bug)
{Mdl = Mat(vector(nb,i,
  my(dl);
  dl = ideallog(nf,RQ[i],bidsmall);
  dl = vectorv(#dl, j, dl[j]%bidcyc[j]);
  dl
))};
*/

/*
makeint(x) = if(type(x)=="t_INT",x,polcoef(x.pol,0));
\\compute discrete logs at large primes with Fp_log_index
install("Fp_log_index","GGGG");
\\g "arith" 1
Lpm = [[p1,m1],[p2,m2],[p3,m3a],[p3,m3b]];
{foreach(Lpm,pm,
  [pi,mi] = pm;
  modpr = nfmodprinit(nf, idealprimedec(nf,pi)[1]);
  g = lift(znprimroot(pi));
  Rmodp = vector(#R,i,makeint(nfmodpr(nf,RQ[i],modpr)));
  print1("computing discrete logs... size m: ", exponent(mi), " size p: ", exponent(pi));
  t = getabstime();
  Ldl = Fp_log_index(Rmodp, g, mi, pi);
  t = getabstime()-t;
  print(" done.");
  printtime(t);
  Mdl = matconcat([Mdl;Ldl]);
  bidcyc = concat(bidcyc,mi);
)};

print1("computing kernel...");
t = getabstime();
[,K] = matsolvemod(Mdl,bidcyc~,0,1);
t = getabstime()-t;
print(" done.");
printtime(t);

print1("LLL reduction...");
t = getabstime();
M = qflll(M*K,3);
t = getabstime()-t;
print(" done.");
printtime(t);
*/

/* check */
/*
H = mathnfmodid(M,h*2^10*(p1-1)*(p2-1)*(p3-1));
hsub = vecprod(vector(#H[,1],i,H[i,i]));
if(hsub != h*4*4*(p1-1)*(p2-1)*(p3-1), error("wrong class number"));

system(strprintf("mv RELSUB-ex2-%d RELSUB-ex2-%d-prev", nb, nb));
write(strprintf("RELSUB-ex2-%d",nb),M);
*/

