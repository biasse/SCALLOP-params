/* relations-1500.gp */

suffix = "1500";

default(parisize,"1G");

d = -4657663250639184493417860902365641674114776564696771899003406891564171427055;
d0 = d \ 9;
factd = factor(d)[,1];

pol = x^2-d0;
nf = nfinit(pol);

addprimes([p | p <- factd, p>10^7]);
a=155853646279738557685783060518469225577;
N=2^126;
fsmall = 3^2 * 32;
p1 = a;         \\ p1,p2,p3
p2 = a-N;       \\ have less than 128 bits
p3 = (a+N)\3;   \\ we can use pari
p4 = a^2-3*N^2; \\ p4 has 250 bits, use cado-nfs
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

system(strprintf("touch STRUCTURE-%s",suffix));
data = readvec(strprintf("STRUCTURE-%s",suffix));
{if(#data,
  print("loading already computed relations for maximal order.");
  h = data[1];
  cyc = data[2];
,\\else
  print("computing relations for maximal order.");
  res = quadclassunit(d0);
  /*
  on paridev:
  cpu time = 3h, 8min, 13,189 ms, real time = 3h, 7min, 27,140 ms.
  */
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


system(strprintf("touch LLLREL-%s",suffix));
Mlll = read(strprintf("LLLREL-%s",suffix));
{if(type(Mlll)=="t_INT",
  print("LLL reducing relations for maximal order.");
  H = mathnfmodid(M,2*h);
  setdebug("qflll",4);
  U = qflll(H,1); \\no flatter
  Mlll = H*U; 
  /* 50 min? TODO check */
  setdebug("qflll",0);
  write(strprintf("LLLREL-%s",suffix), Mlll);
,\\else
  print("loading LLL reduced relations for maximal order.");
)};
\\here we could check that snf(Mlll) is consistent, but this is a bit slow.


nb = 75; \\number of primes we want to keep
forbidden = Set(concat([2,3],factd~));
deleteprimes = List();
keepprimes = List();
{for(i=1,#FB,
  my(p = FB[i]);
  if(setsearch(forbidden,p),
    listput(~deleteprimes,i);
  ,/* else if */#keepprimes<nb,
    listput(~keepprimes,i);
  ,\\else
    listput(~deleteprimes,i);
  );
)};
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
\\could also store/load subFB
if(matsnf(Msub,4) != cyc, error("incompatible Msub SNF!"));


/*
  convention:
  D = -d
  chosen qf of norm p:
  p*x^2 + b*x*y + c*y^2
  with
  b = +- sqrt(D mod p) in [1,p-1]
  s.t. b=D mod 2 (here odd)
  Not sure about even primes (when 2 splits).
*/
dec2 = idealprimedec(nf,2);
orderdec2 = [dec2[2],dec2[1]];
orderdec(p) =
{
  my(dec,b);
  if(p==2,return(orderdec2));
  dec = idealprimedec(nf,p);
  if(#dec==1,return([dec[1],dec[1]]));
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
  id = idealfactorback(nf, vector(nb,i,orderdec(subFB[i])[if(sign(C[i])>0,1,2)]), abs(C));
  Gid = id~*G*id;
  T = qflllgram(Gid);
  res = id*T[,1];
  if(nfeltnorm(nf,res) != factorback(subFB,abs(C)),
    error("wrong norm in reconstruct: ", res)
  );
  res
};

system(strprintf("touch ALGREL-%s-%s",suffix,nb));
R = readvec(strprintf("ALGREL-%s-%s",suffix,nb));
{if(#R==0,
  R = List();
  print("computing algebraic form of the relations.");
  my(L = Vec(Msub));
  print1("reconstruct:");
  t = getabstime();
  for(i=1,nb,
    if(i%20==0,print1(" ",round(100*i/nb),"%"));
    my(C = L[i], g);
    g = reconstruct(C);
    listput(~R,g);
    write(strprintf("ALGREL-%s-%s",suffix,nb), g);
  );
  R = Vec(R);
  t = getabstime()-t;
  print(" done.");
  printtime(t);
,\\else
  print("loading algebraic form of the relations.");
)};

RQ = [[r,1;cc*r,-1] | r <- R];

\\compute discrete logs at small primes with idealstar with cycmod
print("computing discrete logarithms for small primes.");
bidsmall = idealstar(nf,f,1,Msmall);
bidcyc = [gcd(d,Msmall) | d <- bidsmall.cyc]; \\should be done by idealstar (bug)
{Mdl = Mat(vector(nb,i,
  my(dl);
  dl = ideallog(nf,RQ[i],bidsmall);
  dl = vectorv(#dl, j, dl[j]%bidcyc[j]);
  dl
))};

makeint(x) = if(type(x)=="t_INT",x,polcoef(x.pol,0));


install("Fp_log_index","GGGG");
system(strprintf("touch MEDDL-%s-%s",suffix,nb));
Ldl = readvec(strprintf("MEDDL-%s-%s",suffix,nb));
\\Lpm = [[p1,m1a]]; \\for testing
Lpm = [[p1,m1a],[p1,m1b],[p2,m2],[p3,m3a],[p3,m3b]];

print("loading discrete logarithms for medium primes: ", #Ldl, "/", #Lpm);
{for(i=1,#Ldl,
  Mdl = matconcat([Mdl;vdl]);
  bidcyc = concat(bidcyc,Lpm[i][2]);
)};
{if(#Ldl<#Lpm,
  print("computing remaining discrete logarithms for medium primes: ", #Lpm-#Ldl, "/", #Lpm);
  /* ~ 25 minutes */
  setdebug("arith",1);
  for(i=#Ldl+1,#Lpm,
    my([pi,mi] = Lpm[i],modpr,g,Rmodp,t,vdl);
    modpr = nfmodprinit(nf, idealprimedec(nf,pi)[1]);
    g = lift(znprimroot(pi));
    Rmodp = vector(#R,i,makeint(nfmodpr(nf,RQ[i],modpr)));
    print1("computing discrete logs... size m: ", exponent(mi), " size p: ", exponent(pi), " ");
    t = getabstime();
    vdl = Fp_log_index(Rmodp, g, mi, pi);
    t = getabstime()-t;
    print(" done.");
    printtime(t);
    Mdl = matconcat([Mdl;vdl]);
    bidcyc = concat(bidcyc,mi);
    write(strprintf("MEDDL-%s-%s",suffix,nb), vdl);
  );
  setdebug("arith",0);
)};

\\compute discrete logs at large primes with cado-nfs
\\TODO DLP with cado (p4)

/*

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

\\TODO check relations are in suborder by reconstruction

H = mathnfmodid(M,h*2^10*(p1-1)*(p2-1)*(p3-1));
hsub = vecprod(vector(#H[,1],i,H[i,i]));
if(hsub != h*4*4*(p1-1)*(p2-1)*(p3-1), error("wrong class number"));

system(strprintf("mv RELSUB-ex2-%d RELSUB-ex2-%d-prev", nb, nb));
write(strprintf("RELSUB-ex2-%d",nb),M);
*/

