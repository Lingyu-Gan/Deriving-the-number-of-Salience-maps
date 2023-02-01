% this script is to plot figure 5

% Fig 5a % uncorrect multi-item
run IdealDetectorWithPerfAttenFilter.m
run IdealDetectorWithSubjectAttenFilter.m
% Fig 5b % corrrected with motor error
run getMotorError.m
run CorrectForMotor_PerfectFilter.m
run CorrectForMotor_subjectAttenFiter.m
