add_custom_target(MaCh3DUNEApps)

foreach(app
	EventRates
	Fit
	# Projections
	SigmaVariation
	LikelihoodScan
	# Variations
	SplineMaker
	# SigmaVarTutorial
	MakeCovarianceMatrix
      )

    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} SamplesDUNE MaCh3DUNE::All Minuit Eigen3::Eigen)

    add_dependencies(MaCh3DUNEApps ${app})

    install(TARGETS ${app} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach(app)

# add_subdirectory(AbiP/Plotting)
