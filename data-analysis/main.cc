/******************************************************************************
 *
 * This is an example main()
 *
 *****************************************************************************/

#include "Wasa.hh"
#include "SorterConfig.hh"

int main(int argc, char** argv) {
  // read command line arguments, mandatory for proper operation
  gSorterConfig->ReadCmdLine(argc,argv);

  // register Sorter ID
  Wasa::Initialize("eventselection","","RootSorter.log");

  // register analysis module
  gWasa->AddAnalysis("eventselection","Raw");

  gWasa->Run();

  delete gWasa;
}



