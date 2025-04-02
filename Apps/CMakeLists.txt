add_custom_target(MaCh3DUNEApps)

foreach(app
	EventRates
	Fit
	Projections
	SigmaVariation
	LikelihoodScan
  Variations
      )
      
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} MaCh3DUNE::All Minuit)

    add_dependencies(MaCh3DUNEApps ${app})

    install(TARGETS ${app} DESTINATION bin)
endforeach(app)
