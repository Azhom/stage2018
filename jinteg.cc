/* ----   
   Exemple d'integration numerique avec SOPHYA/NTools 
   pour Jeremy  -   Fichier jinteg.h 
   execution par runcxx :
   runcxx -inc jinteg.h -f jinteg.cc
*/

// On instancie nos classe-fonction
myFuncPoly f(2., 3.);
cout << " ---- test de myFuncPoly f(3., 3.); " << endl;
for(double x=-2.; x<12.; x+=2.) {
  cout<<" x="<<x<<" f(x)="<<f(x)<<endl;
 }
myFuncGauss g(5., 1.);
cout << " ---- test de myFuncGauss g(5., 1.) " << endl;
for(double x=-2.; x<12.; x+=2.) {
  cout<<" x="<<x<<" g(x)="<<g(x)<<endl;
 }

//  On calcule l'integrale par la methode des Trapezes
//  (Note: Il faudra ameliorer l'interface et fonctionalites de la classe) 
double xmin=-10.; double xmax=50.;  // limites d'integration
TrpzInteg trpinteg;
// Nb de pas de d'integration
int npt=1000;
trpinteg.NStep(npt);   
trpinteg.SetFunc(f);
cout << " ----Integrale_Trapeze[f, xmin, xmax]="<<trpinteg.ValueBetween(xmin, xmax)<<endl;
trpinteg.SetFunc(g);
cout << " ----Integrale_Trapeze[g, xmin, xmax]="<<trpinteg.ValueBetween(xmin, xmax)<<endl;

//  On calcule l'integrale par la methode des Gauss-Legendre
//  (Note: Il faudra ameliorer l'interface et fonctionalites de la classe) 
double x1=0.; double x2=10.;  // limites d'integration
GLInteg glinteg;
// On definit l'ordre (degre) des polynoms 
int order=200;
glinteg.SetOrder(order);  // Nb de pas de d'integration
glinteg.SetFunc(f);
cout << " ----Integrale_GaussLegendre[f, xmin, xmax]="<<glinteg.ValueBetween(xmin, xmax)<<endl;
glinteg.SetFunc(g);
cout << " ----Integrale_GaussLegendre[g, xmin, xmax]="<<glinteg.ValueBetween(xmin, xmax)<<endl;
