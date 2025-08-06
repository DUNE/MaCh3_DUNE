add_custom_target(MaCh3DUNEApps)

foreach(app
	EventRates
	Fit
	Projections
	SigmaVariation
	LikelihoodScan
	Variations
	EventRates_New
	Fit_Abi
	Fit_fixingparameters
	Fit_fixingparameters_beluga
	Fit_fixingparameters_beluga_withthrows
	Abi_Posteriorpredicitive
	Abi_Posteriorpredictive_fardetector
	Postfit_xsec_diagnostics
      )
      
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} MaCh3DUNE::All Minuit)

    add_dependencies(MaCh3DUNEApps ${app})

    install(TARGETS ${app} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach(app)
