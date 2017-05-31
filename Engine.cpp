

// RooBarb
#include "Utils.h"
#include "Logger.h"
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

	cout << "test: "  << quote( "hello" ) << endl;
	Logger::setGlobalLogLevel( "none" );
	// XmlConfig cfg( argv[1] );

	TaskFactory::registerTaskRunner<PicoDstSkimmer>( "PicoDstSkimmer" );
	TaskFactory::registerTaskRunner<FemtoDstWriter>( "FemtoDstWriter" );

	// XmlConfig cfg( argv[1] );
	TaskEngine engine( argc, argv );

	return 0;
}
