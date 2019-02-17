SET(${PROJECT_NAME}_HELPERS_HEADERS
  definitions.h
  details.h
  linear_problem.h
  quadratic_cost.h
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HELPERS_HEADERS}
  DESTINATION include/hpp/spline/optimization
  )
