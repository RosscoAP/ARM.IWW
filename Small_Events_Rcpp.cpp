#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

/*
C++ function to add small events to single stations
*/
// [[Rcpp::export]]
IntegerVector insert_small_events_station(
	IntegerVector Pcp, // vector of Pcp in 1/100th mm
	IntegerVector Rain, // logical vector whether it's raining or not
  IntegerVector EventIdxs, // event indexes for the current season
	IntegerVector EventStart, // starting timestep for each dry spell
	IntegerVector EventL, // length of dry spells
	IntegerMatrix SmallEvents // columns: WSD (0), WSI (1)
) {
  IntegerVector StartIdx(SmallEvents.nrow()); // create vectors to store starting position for each small event.
  NumericVector RelDistance(SmallEvents.nrow()); // and the relative starting position from DSD start
  IntegerVector DSDIdx(SmallEvents.nrow()); // and the relative starting position from DSD start
  bool success = 0; // just a switch to designate whether a small event was successfully placed or not
  bool failed = 0; // a similar switch used for the same purpose

  // loop around small events and place!
    for (int i = 0; i < SmallEvents.nrow(); i++) {
      success = 0; // reset switch for looping
			//Rcout <<  "i: " << i << " of " << SmallEvents.nrow() << "\n";
      while (!success) {
        failed = 0; // reset switch
        // randomly select a dry spell from the season
        DSDIdx(i) = sample(EventIdxs, 1)(0);
				//Rcout <<  "DSDIdx: " << DSDIdx(i) << "\n";
        // check if DSD is large enough (WSD + 2)
        if (EventL(DSDIdx(i)) <  SmallEvents(i,0) + 2  ) {continue;} // start again
        // select a starting position within the DSD
        RelDistance(i) = R::runif(0, 1);
        StartIdx(i) = round( EventStart(DSDIdx(i)) + 1 + RelDistance(i) * (EventL(DSDIdx(i)) - 2 -  SmallEvents(i,0) ) );
        //Rcout <<  "i: " << i << " of " << SmallEvents.nrow() << " StartIdx: " << StartIdx(i) << " RelDistance: " << RelDistance(i) << "\n";
				//Rprintf("i: %i StartIdx: %i \n", i, StartIdx(i) );

        // check if already occupied by a small event
        for(int j = -1; j < SmallEvents(i,0) + 1; j++) { // loop around timesteps of small event (including timestep before and after)
          if (Rain(StartIdx(i)+j)==1) { failed=1; break; }
        }
        if (failed) { continue; } // try again

        // should be good now to place the small event
        for(int j = 0; j < SmallEvents(i,0); j++) { // loop around timesteps of small event
          Pcp(StartIdx(i) + j ) += SmallEvents(i,1);
          Rain(StartIdx(i) + j) = 1;
        }
        // all good, continue to next small event
        success = 1; // to stop looping
      }
    }
  // only return pcp vector
  return Pcp; // for some reason it does this strange thing where the Pcp is updated automatically in the workspace. You can't save it to a different variable
}
