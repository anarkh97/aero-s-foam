ADD_LIBRARY( pita 
		  LinearizedTimeIntegrator.C
		  NLDistrTimeDecompSolver.C
		  NLDynamTimeIntegrator.C
		  NLTimeSlice.C
		  NiceTimer.C
		  PitaNonLinDynam.C
		  TimeSliceMapping.C
                  ../Driver.d/DomainPita.C
)
