// Für Ein-Ausgabe
#include <iostream>
#include <fstream>
#include <string>
// Definitionen verschiedener Exception Klassen. Hier gebraucht für runtime_error
#include <stdexcept>
// Boost program_options macht das Handhaben von Optionen leichter.
#include <boost/program_options.hpp>
// c++ version von printf.
#include <boost/format.hpp>
// ublas bietet verschiedene Klassen für Lineare Algebra. Siehe unten.
#include <boost/numeric/ublas/vector.hpp>


namespace ublas = boost::numeric::ublas;

using std::string;
using std::cout;
using std::endl;


namespace PredatorPreyPdeIntegration
{
/*
Variablen für die Modellparameter sind im Namespace "P" deklariert. Die
Variablennamen sind absichtlich kurz gehalten und analog zur Modelldefinition.

Zur Unterscheidung von anderen Variablen dient der Namespace. Ganz allgemein
bezieht man sich auf "Etwas" (kann Variable, Funktion oder Klasse sein) in
einem Namespace mithilfe des "::" Operators. Beispiele:
  std::cout << "bla bla" << std::endl.
oder
  std::vector<int> intvector;

Man kann auch durch "using namespace Name", also z.b "using namespace std" alle
Objekte aus einem Namespace importieren so, daß man Name:: nicht mehr braucht.

Der gesammte restliche Code ist ebenfalls in einem Namespace, nämlich
PredatorPreyPdeIntegration. Zu einem bestimmten Namen sucht der Compiler das
passende Objekt stets zuerst im "eigenen" Namespace. Das vermeidet
Zweideutigkeiten. Ohne den PredatorPreyPdeIntegration Namespace gäb es einen
Compilerfehler, weil in einem der inkludierten Headers bereits der Name "time"
für etwas vergeben wurde, und dieser Name hier auch für eine Variable verwendet wird.
*/
namespace P // Parameters
{
double d1 = 1.,
       d2 = 1.,
       r = 1.,
       w = 1.,
       p = 1.,
       k = 1.,
       q = 1.,
       s = 1.;
int init = 0;
}
/*
Kurz: Zugriff auf Modellparameter q durch P::q. z.b. a = b * P::q + ...
*/


void print_description(std::ostream &os)
{
  os << "the spatial predator-prey-model:" << endl;
  os << "  x' =  d1 laplace x + r*x*(1-x/w) - p*y*h(k*x)" << endl;
  os << "  y' =  d2 laplace y + q*h(k*x) - s*y" << endl;
  os << "  h = eta -> eta/(1+eta)" << endl;
  os << "  x = prey, y = predators" << endl;
}

/*
PRED und PREY werden definiert als Namen für 0 bzw 1.
Sie dienen nur der Klarheit. Man könnte auch 0 oder 1 schreiben.
NUM_SPECIES steht einfach für die Anzahl der Spezies.
*/
enum {
  PRED = 0,
  PREY = 1,
  NUM_SPECIES = 2
};

const char* species_names[NUM_SPECIES] = {
  "predators",
  "prey"
};

/*
Für jeden Gitterpunkt muss der Wert der Spezieskonzentration gespeichert werden.
Dafür nehmen wir die vector Klasse aus der Ublas Bibliothek. Intern nutzt sie
ein dynamisches Array, wie z.b. double *daten = new double[number_of_sites],
erleichtert aber das Speichermanagement und unterstützt arithmetische Operationen.

typedef ist dazu da einem Datentyp einen zusätzlichen Namen zu geben. Wird sehr
häufig genutzt um kürzere Schreibweisen zu ermöglichen.
*/
typedef ublas::vector<double> Vector;

/*
 SpeciesArray wird nun den Typ "statisches Array aus NUM_SPEZIES Vektoren" bezeichnen.
 Im Gegensatz zum ublas::vector steht die Länge dieses Arrays zur Zeit der
 Kompilation fest.
 
 Wir können SpeciesArrays als Argumente an Funktionen übergeben (Bsp: get_initial_value 
 weiter unten). Um an die Daten zu kommen können wir z.b.
 double predators_an_stelle_0 = solution[PRED][0]
 schreiben, wobei solution vom Typ SpeciesArray ist und PRED vorher als 0 definiert wurde.
 Wir können aber auch per Schleife über die Spezies iterieren.
*/
typedef Vector SpeciesArray[NUM_SPECIES];


// Allgemeine Programm Optionen
double max_time = 100.,
       dt = 0.1,
       output_intervall = 0.1,
       spacing = 1.;
int    grid_size = 50;   
string output_filename = "pred-prey-pde-testrun.dat";
// Aktueller Wert zur Zeit t
SpeciesArray solution;
// "Buchhalte"-daten
int num_sites;
int output_number;
double time, next_output_time;



inline double gauss_packet(double x, double y, double w)
{
  // Berechnet den Wert eines Gausspackets mit der Breite w.
  return std::exp(-(x*x+y*y)/(w*w));
}


void fill_with_initial_value(SpeciesArray &v)
{
  int c = 0; // c ist die Nummer des aktuellen Gitterplatzes.
  // Sie wird unten in der zweiten "for" Schleife zusammen mit "x" jede Iteration um 1 erhöht.
  // Die Nummerierung der Gitterplätze muss natürlich über das ganze Programm konsistent sein.
  for (int y=0; y<grid_size; ++y)
    for (int x=0; x<grid_size; ++x, ++c)
    {
      switch(P::init)
      {
        case 0:
          // Konfiguration 0 ist jeweils eine "Rampen"-Funktion, einmal in x und einmal in y Richtung.
          v[PRED][c] = 0.02 * double(grid_size - y) / grid_size + 0.31;
          v[PREY][c] = 0.02 * double(x) / grid_size + 0.19;
          break;
        case 1:
        { // Scope aufgemacht damit man lokale Variablen deklarieren kann.
          // Konfiguration 1 ist jeweilt ein Gausspacket.
          // Es folgen die Koordinaten der Packetmittelpunkte.
          double w = grid_size * 0.1;
          double pred_x = 0.5*(grid_size - w);
          double pred_y = 0.5*(grid_size);
          double prey_x = pred_x + w;
          double prey_y = pred_y;
          // Und nun die Funktion
          v[PRED][c] = 0.5 * gauss_packet(x-pred_x, y-pred_y, w);
          v[PREY][c] = 0.5 * gauss_packet(x-prey_x, y-prey_y, w);
        }
        break;
      }
    }
}


void add_discrete_laplacian(Vector &dst, const Vector &v, double prefactor)
{
  // Implement me! :)

  /*
  Hinweise:
  "dst" und "v" werden per Referenz übergeben (bedingt durch das &). Das bedeutet,
  daß keine Kopie gemacht wird. Es wird lediglich ein Verweis ähnlich wie ein Zeiger
  auf die übergeben Daten verwendet. Diese können also innerhalb der Funktion verändert werden.

  "v" ist zusätzlich als const deklariert. Dadurch ist nur eine Untermenge der normalen
  Operationen zugänglich, im Normalfall solche, die das Objekt unverändert lassen.
  
  Schauen Sie wie im Rest des Programms über die Gitterplätze iteriert wird.
  
  Am einfachsten ist es wohl no-flux Randbedingungen mit Hilfe von Bedingungen wie
  "if (x > 0) ...;  if (x < grid_size-1) ...;" zu implementieren. So kann man diffusive
  Flüsse zwischen benachbarten Plätzen aufaddieren bzw. an Randpunkten entsprechend weglassen.
  */
}


void integrate()
{
  // Temporärer Speicher. 
  SpeciesArray f = {
    // Initialisiert je einen ublas::vector mit der richtigen Länge und Nullen.
    ublas::scalar_vector<double>(num_sites, 0.),
    ublas::scalar_vector<double>(num_sites, 0.)
  };
  /*
  Lösungsentwurf
  (*) Mittels add_discrete_laplacian diffusive Beiträge in f speichern.
  (*) Lokale Beiträge von Quellen-und Senktenterme auf f addieren.
  (*) Werte in "solution" updaten
  */

  // Implement me! :)

  time += dt;
}


void init()
{
  num_sites = grid_size * grid_size;
  time = next_output_time = 0.;
  output_number = 0;

  double stability_check_number = 4*dt/(spacing*spacing)*std::max(P::d1, P::d2);
  if (stability_check_number > 1.)
    cout << "Warning: Stability condition violated. 4*D*dt*spacing^2 = "
         << stability_check_number << " must be < 1" << endl;

  solution[PRED].resize(num_sites);
  solution[PREY].resize(num_sites);
  fill_with_initial_value(solution);
}


void output_results()
{
  /*
  * Hier geben wir die aktuelle Lösung im Gnuplot-tauglichen Format aus:
  * Je Zeile: x y spezies1 spezies2
  * x variiert am schnellsten. Nach jeder Änderung von y kommt eine Leerzeile.
  * Zeilen, die mit # anfangen sind Kommentare.
  */
  string fn = (boost::format("%s-%03i.dat")
                  % output_filename
                  % output_number).str();

  std::ofstream f(fn.c_str());
  if (!f.good())
    throw std::runtime_error("unable to open output file: "+fn);

  // Wir rechnen noch die Mittelwerte der Konzentrationen aus.
  double avg_values[NUM_SPECIES];
  for (int s=0; s<NUM_SPECIES; ++s)
  {
    avg_values[s] = 0.;
    for (int i=0; i<num_sites; ++i)
      avg_values[s] += solution[s][i];
  }
  f << "#time=" << time << endl;
  for (int s=0; s<NUM_SPECIES; ++s)
  {
    avg_values[s] /= num_sites;
    f << "#avg-" << species_names[s] << "=" << avg_values[s] << endl;
  }

  // Anzahl der Gitterplätze
  f << boost::format("#grid-size=%i") % grid_size << endl;

  // Hier werden die Gitterwerte geschrieben.
  int i = 0;
  for (int y=0; y<grid_size; ++y)
  {
    for (int x=0; x<grid_size; ++x, ++i)
    {
      f << x << " " << y;
      for (int s=0; s<NUM_SPECIES; ++s)
      {
        f << " " << solution[s][i];
      }
      f << endl;
    }
    f << endl;
  }

  // Check ob noch alle Werte endlich sind. Sonst Programmabbruch.
  string status = (boost::format("@time %f") % time).str();
  for (int s=0; s<NUM_SPECIES; ++s)
  {
    if (!std::isfinite(avg_values[s]))
      throw std::runtime_error("non-finite solution values");
    status += (boost::format(", <%s>=%f") % species_names[s] % avg_values[s]).str();
  }
  cout << status << endl;
  ++output_number;
}


void run()
{
  init();
  while (true)
  {
    if (time >= next_output_time - 1.e-3 * dt)
    {
      output_results();
      next_output_time += output_intervall;
    }

    if (time > max_time) break; // fertig

    integrate();
  }
}



void handle_program_args(int argc, char **argv)
{
  /* boost program_options basierter "Kochplattencode" für Programmargumenthandling
   * und hübsche Usage Anweisungen (--help).
   * 
   * Parametervariablen im P Namespace werden hier gesetzt!
   */
  
  namespace po = boost::program_options;

  // Makro für Program Optionen in die options_description Liste zu machen.
  #define ADD_PARAM(NAME, P, MSG) (NAME, po::value(&P)->default_value(P), MSG)
  po::options_description desc("Options");
  desc.add_options()
    ("help", "generate help message")
    ADD_PARAM("mp-d1", P::d1, "prey diffusion coefficient")
    ADD_PARAM("mp-d2", P::d2, "predator diffusion coefficient")
    ADD_PARAM("mp-r", P::r, "prey proliferation coefficient")
    ADD_PARAM("mp-w", P::w, "maximum prey density")
    ADD_PARAM("mp-p", P::p, "prey death coefficient")
    ADD_PARAM("mp-k", P::k, "prey scaling coefficient")
    ADD_PARAM("mp-q", P::q, "predator proliferation coefficient")
    ADD_PARAM("mp-s", P::s, "predator death coefficient")
    ADD_PARAM("mp-initcase", P::init, "select an initial species distribution. Values are 0 = ramp and 1 = gauss packet.")
    ADD_PARAM("max-time", max_time, "simulation duration")
    ADD_PARAM("time-step", dt, "time step")
    ADD_PARAM("grid-size", grid_size, "number of grid points along each dimension")
    ADD_PARAM("grid-spacing", spacing, "grid spacing")
    ADD_PARAM("out-fn", output_filename, "output filename")
    ADD_PARAM("out-intervall", output_intervall, "output intervall" )
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help"))
  {
    cout << desc << endl;
    cout << "The PDE is";
    print_description(cout);
    exit(1);
  }
  cout << "Integrating "; print_description(cout);
  cout << "Model Parameters:" << endl;
  #define OUT_PARAM1(PARAM) cout << #PARAM" = " << P::PARAM << endl
  OUT_PARAM1(d1);
  OUT_PARAM1(d2);
  OUT_PARAM1(r);
  OUT_PARAM1(w);
  OUT_PARAM1(p);
  OUT_PARAM1(k);
  OUT_PARAM1(q);
  OUT_PARAM1(s);
  cout << boost::format("Lattice: %ix%i, spacing: %f") % grid_size % grid_size % spacing << endl;
  cout << boost::format("max-time: %f, dt = %f, output-interval = %f") % max_time % dt % output_intervall << endl;
  cout << "output-filename: " << output_filename << endl;
}

}


int main(int argc, char **argv)
{
  PredatorPreyPdeIntegration::handle_program_args(argc, argv);
  PredatorPreyPdeIntegration::run();
  return 1;
}
