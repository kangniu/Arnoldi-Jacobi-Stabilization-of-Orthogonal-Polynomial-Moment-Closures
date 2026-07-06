function Uclose = closure_raw_moment(U,closure)
%CLOSURE_RAW_MOMENT Compute U_{n+1} from known moments U_0,...,U_n.

switch upper(closure)
    case 'B4'
        % Known U0,...,U3.  Closing normalized central moment:
        % W4 = 2 W3^2 + 3.
        [rho,u,theta] = raw_to_primitive(U);
        W = raw_to_normalized(U,3);
        Wfull = zeros(5,1);
        Wfull(1) = 1;
        Wfull(2) = 0;
        Wfull(3) = 1;
        Wfull(4) = W(4);
        Wfull(5) = 2*W(4)^2 + 3;
        Uclose = central_normalized_to_raw(rho,u,theta,Wfull,4);

    case 'HYQMOM5'
        % Known U0,...,U4. Closing W5:
        % W5 = 1/2 W3(5 W4 - 3 W3^2 - 1).
        [rho,u,theta] = raw_to_primitive(U);
        W = raw_to_normalized(U,4);
        Wfull = zeros(6,1);
        Wfull(1) = 1; Wfull(2) = 0; Wfull(3) = 1;
        Wfull(4) = W(4);
        Wfull(5) = W(5);
        Wfull(6) = 0.5*W(4)*(5*W(5) - 3*W(4)^2 - 1);
        Uclose = central_normalized_to_raw(rho,u,theta,Wfull,5);

    case 'B6'
        % Known U0,...,U5. Closing W6 formula from Morin--McDonald B6.
        [rho,u,theta] = raw_to_primitive(U);
        W = raw_to_normalized(U,5);
        W3 = W(4);
        W4 = W(5);
        W5 = W(6);

        % B6 closing flux derived directly from the closure coefficient
        %
        %   b3 = (b1+b2) + 1/2[(a2-a0)^2 + (a2-a1)^2].
        %
        % This form satisfies the equilibrium check:
        %   (W3,W4,W5)=(0,3,0) -> W6=15.
        %
        % It is algebraically equivalent to solving the recurrence condition
        % for W6.  This is safer than relying on visually parsed PDF signs.
        g = W3^2 - W4 + 1;
        if abs(g) < 1e-10
            g = sign_nonzero(g)*1e-10;
        end

        N6 = 3*W3^6 - 10*W3^4*W4 - 2*W3^4 + ...
             6*W3^3*W5 + 7*W3^2*W4^2 + 8*W3^2*W4 - ...
             W3^2 - 14*W3*W4*W5 - 2*W3*W5 + ...
             4*W4^3 - 6*W4^2 + 2*W4 + 4*W5^2;

        W6 = -0.5*N6/g;

        Wfull = zeros(7,1);
        Wfull(1) = 1; Wfull(2) = 0; Wfull(3) = 1;
        Wfull(4) = W3; Wfull(5) = W4; Wfull(6) = W5; Wfull(7) = W6;
        Uclose = central_normalized_to_raw(rho,u,theta,Wfull,6);

    case 'B8'
        % Known U0,...,U7.  The B8 closing moment is recovered from the
        % hierarchy coefficient b4, avoiding a hand-written W8 formula.
        [rho,u,theta] = raw_to_primitive(U);
        W = raw_to_normalized(U,7);
        W3 = W(4);
        W4 = W(5);
        W5 = W(6);
        W6 = W(7);
        W7 = W(8);

        W8 = b8_closing_moment_normalized(W3,W4,W5,W6,W7);

        Wfull = zeros(9,1);
        Wfull(1) = 1; Wfull(2) = 0; Wfull(3) = 1;
        Wfull(4) = W3; Wfull(5) = W4; Wfull(6) = W5;
        Wfull(7) = W6; Wfull(8) = W7; Wfull(9) = W8;
        Uclose = central_normalized_to_raw(rho,u,theta,Wfull,8);

    otherwise
        error('Unknown closure %s',closure);
end

end

function s = sign_nonzero(x)
if x >= 0
    s = 1;
else
    s = -1;
end
end
