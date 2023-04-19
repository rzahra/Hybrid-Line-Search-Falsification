function res = solve_hybrid(this)
% Solve_quasi_random works with quasi-random sampling   
% 

%seed = this.solver_options.rand_seed;

seed = 
num_samples =  this.solver_options.num_rand_samples;

if ~strcmp(this.display,'off')
    fprintf('\n\n++++++++++++++++++++++++++++++++++++++++++++\nRunning %g corners-random samples with seed %g\n', num_samples, seed);
    this.display_status_header();
end

BrR = this.BrSet.copy();
BrR.ResetParamSet();
BrR.SetParamRanges(this.params, [this.lb this.ub])
rng(seed);
BrR.SampleDomain(this.params,num_samples);
X1 = BrR.GetParam(this.params);







BrC = this.BrSet.copy();
BrC.ResetParamSet();
BrC.SetParamRanges(this.params, [this.lb this.ub])
num_corners =  this.solver_options.num_corners;
BrC.CornerSample(num_corners);
X2 = BrC.GetParam(this.params);



starting = "corners";

down_corners = "No";

corners_counter = 1;
random_counter = 1;

for i = 1:1000
    
    if strcmp (starting, "corners")
        
        res(:,i) = this.FevalInit(X2(:,corners_counter));
        
        starting = "randoms";
        
        
        if corners_counter == length (X2(1,:))
            down_corners = "Yes";
            
        else
            corners_counter = corners_counter + 1;
        end
        
        
        
    elseif strcmp (starting, "randoms")
        
        
        res(:,i) = this.FevalInit(X1(:,random_counter));
        
        if strcmp (down_corners, "No")
            
            starting = "corners";
            
        else
            starting = "randoms";
        end
        
        
        random_counter = random_counter + 1;
        
    end
    
    
end


%res = this.FevalInit(X0);
%this.add_res(res);



end