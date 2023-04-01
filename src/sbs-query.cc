#include "config.h"

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"

#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "observatory.h"
#include "sbsearch.h"
#include "test_db.h"
#include "util.h"

using sbsearch::Indexer;
using sbsearch::Logger;
using sbsearch::Observation;
using sbsearch::Observatories;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
