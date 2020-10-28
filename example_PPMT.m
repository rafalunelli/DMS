
load('datasets/HardData_ReferenceModel_size100_range20.mat');

I_ref = size(reference_models,2);
J_ref = size(reference_models,3);
reference_variables = [reshape(reference_models(1,:,:),I_ref*J_ref,1) reshape(reference_models(2,:,:),I_ref*J_ref,1) reshape(reference_models(3,:,:),I_ref*J_ref,1) reshape(reference_models(4,:,:),I_ref*J_ref,1) reshape(reference_models(5,:,:),I_ref*J_ref,1) reshape(reference_models(6,:,:),I_ref*J_ref,1) ];

I = I_ref;
J = J_ref;
range = 18;
n_simulations = 2;

n_cond_points = 20;
cond_value_ = cond_value(1:n_cond_points ,:);
cond_pos_ = cond_pos(1:n_cond_points ,:);

% Uncondicional PPMT
simulations_all_ppmt_unconditional = run_PPMT(I,J,I_ref,J_ref, [range range], reference_variables, [], [], n_simulations);

% Condicional PPMT
simulations_all_ppmt_conditional = run_PPMT(I,J,I_ref,J_ref, [range range], reference_variables, cond_pos_, cond_value_, n_simulations);


simulation_dms = simulations_all_ppmt_conditional{1};

generate_2D(reference_models,cond_pos_)
set_caxis_6variate
generate_2D(simulation_dms,cond_pos_)
set_caxis_6variate

generate_histograms(reshape(reference_models,6,I*J)')
generate_histograms(reshape(simulation_dms,6,I*J)')

%generate_isotropic_variograms(0.5,reference_models,simulation_dms,reshape(simulation_ppmt',6,I,I))

load('datasets/reference_data.mat')
reference = [z1_analytic z2_analytic, z3_analytic, z4_analytic, z5_analytic, z6_analytic];
num_of_bins = 50;
aux_dms = reshape(simulation_dms,6,I*J)';
%aux_dms = [ simulation(:,1) simulation(:,2) simulation(:,3) simulation(:,4) simulation(:,5) simulation(:,6)];
chi2_dms = generate_chi2(reference,aux_dms, num_of_bins,0);


    



