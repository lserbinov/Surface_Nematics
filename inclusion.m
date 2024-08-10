
n_A = 10;
I_A = (1:n_A)'; % I_A: A -> A identity map
x = linspace(0,1,n_A)';  % x: A -> R
u = sin(x*2*pi);
isB = abs(u)<0.5;
j_B = find(isB); % j_B : B -> A  inclusion map
n_B = length(j_B);
I_B = (1:n_B)'; % : B-> B identity map
j_B_inv = zeros(n_B,1);
j_B_inv(j_B) = I_B;