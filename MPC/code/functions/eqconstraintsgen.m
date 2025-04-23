function eqconstraints = eqconstraintsgen(LTIe, dim, dime, dtilde)
eqconstraints.A = [eye(dim.nx) - LTIe.A(1:dim.nx, 1:dim.nx), -LTIe.B(1:dim.nx, :); 
                                    LTIe.C(:, 1:dim.nx), zeros(dime.ny, dim.nu) ];
eqconstraints.b = [LTIe.A(1:dim.nx, (dim.nx+1):end)*dtilde; 
                     LTIe.yref - LTIe.C(:, (dim.nx+1):end)*dtilde];
end