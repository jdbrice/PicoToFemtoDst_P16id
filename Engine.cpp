

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>


#include "PicoDstSkimmer/PicoDstSkimmer.h"
#include "FemtoDstWriter/FemtoDstWriter.h"

#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);


	TaskFactory::registerTaskRunner<PicoDstSkimmer>( "PicoDstSkimmer" );
	TaskFactory::registerTaskRunner<FemtoDstWriter>( "FemtoDstWriter" );

	TaskEngine engine( argc, argv );

	return 0;
}
