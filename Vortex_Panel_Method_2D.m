%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITTEN BY NIKOLAI SMITH, MAY 2024                                                            % 
% The following program is a vortex panel method solver for the 2D                              %
% incompressible, inviscid potential flow problem.                                              %
% In this solver:                                                                               %
%  - NACA airfoils are discretized into distinct constant-strength vortex                       %
%    panels                                                                                     %
%  - The internal tangential influence coefficients from all                                    %
%    panels at a point directly under the current panel is computed for all                     %
%    panels in the problem                                                                      %
%  - The vortex strengths are selected to satisfy the indirect boundary                         %
%    condition problem (Dirichlet BC)                                                           %
%  - After computing the vortex strengths, the velocity field is computed                       %
%    at each collocation point along the airfoil surface to acquire the                         %
%    pressure coefficients which can be numerically integrated to acquire                       %
%    the lift coefficient of that airfoil                                                       %
%  - Additional supplementary functions are added to create pressure and                        %
%    flow field visualization.                                                                  %
% REFERENCES:                                                                                   %
%   LOW SPEED AERODYNAMICS:                                                                     %
%       Katz J, Plotkin A. Low-Speed Aerodynamics. 2nd ed. Cambridge University Press; 2001.    %
%   JoshTheEngineer:                                                                            %
%       https://www.youtube.com/user/JoshTheEngineer                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all; %#ok
n_panels = 2000;
max_camber = 2.0; % Floating point between 0 and 10
max_camber_loc = 4.0; % Floating point number between 0 and 10
max_percent_thickness = 12.0; % Floating point number between 10 and 30
attack_angle = 10.0; % Floating point number between -14 and 14

coordinates = gen_airfoil_coords( max_camber , max_camber_loc , max_percent_thickness , n_panels );
[ norm_rel_velos, tan_rel_velos ] = gen_relative_velos( coordinates, attack_angle );
[ norm_influence_coeffs , tan_influence_coeffs ] = gen_boundary_influence_coeffs( coordinates );
vortex_strengths = calculate_vortex_strengths( tan_influence_coeffs , tan_rel_velos );
pres_coeffs = pressure_coeffs( tan_influence_coeffs , vortex_strengths , tan_rel_velos );
std_lift = standard_lift( coordinates , pres_coeffs );
kj_lift = kutta_joukowski_lift( coordinates , vortex_strengths );
fprintf("Standard Lift Coefficient: %f\n", std_lift);
fprintf("Kutta-Joukowski Lift Coefficient: %f\n", kj_lift);
plot_pres_coeffs( coordinates , pres_coeffs );
%streamline_plot( coordinates , vortex_strengths , attack_angle );
%pressure_contour_plot( coordinates , vortex_strengths , attack_angle );

function coordinates = gen_airfoil_coords( m , p , t , n_panels )

    m = m / 100;
    p = p / 10;
    t = t / 100;

    n_eval_points = n_panels + 1;
    eval_points = ( 1 + cos( linspace( 0 , 2*pi , n_eval_points ) ) ) / 2;
    x_coords = zeros( [ n_eval_points 1 ] );
    y_coords = zeros( [ n_eval_points 1 ] );
    [ ~ , min_index ] = min( eval_points );

    for i = 1 : n_eval_points

        x = eval_points( i );

        y_t = 5 * t * ( 0.2969 * sqrt( x ) - 0.1260 * x - 0.3516 * x ^ 2 + 0.2843 * x ^ 3 - 0.1015 * x ^ 4 );

        if x < p

            y_c = m / p ^ 2 * ( 2 * p * x - x ^ 2 );
            y_c_deriv = 2 * m / p ^ 2 * ( p - x );

        else

            y_c = m / ( 1 - p ) ^ 2 * ( ( 1 - 2 * p ) + 2 * p * x - x^2 );
            y_c_deriv = 2 * m / ( 1 - p ) ^ 2 * ( p - x );

        end

        theta = atan( y_c_deriv );
        if ( i < min_index + 1 )
            x_coords( i ) = x - y_t * sin( theta );
            y_coords( i ) = y_c + y_t * cos( theta );
        else
            x_coords( i ) = x + y_t * sin( theta );
            y_coords( i ) = y_c - y_t * cos( theta );
        end

    end

    coordinates = zeros( [ n_eval_points 2 ] );
    coordinates( : , 1 ) = x_coords;
    coordinates( : , 2 ) = y_coords;

end

function [ norm_freestream_velos , tan_freestream_velos ] = gen_relative_velos( coordinates , attack_angle )

    attack_angle = attack_angle * pi / 180;
    freestream_velo = [ cos( attack_angle ) ; sin( attack_angle ) ];
    n_panels = length( coordinates( : , 1 ) ) - 1;
    norm_freestream_velos = zeros( [ n_panels 1 ] );
    tan_freestream_velos = zeros( [ n_panels 1 ] );

    for i = 1 : n_panels

        panel = coordinates( i + 1 , : ) - coordinates( i , : );

        panel_norm = [ - panel( 2 ) ; panel( 1 ) ] / norm( panel );
        panel_tan = [ panel( 1 ) ; panel( 2 ) ] / norm( panel );

        norm_freestream_velos( i ) = dot( freestream_velo , panel_norm );
        tan_freestream_velos( i ) = dot( freestream_velo , panel_tan );

    end

end

function world_influence_coeffs = gen_world_influence_coeffs( coordinates , x_p , y_p )
    n_panels = length( coordinates( : , 1 ) ) - 1;
    world_influence_coeffs = zeros( [ 2 n_panels ] );
    for i = 1 : n_panels
        if ( ( coordinates( i , : ) + coordinates( i + 1 , : ) ) / 2 ~= [ x_p y_p ] )
                x_k = coordinates( i , 1 );
                y_k = coordinates( i , 2 );
                panel = coordinates( i + 1 , : ) - coordinates( i , : );
                phi = atan2( panel( 2 ) , panel( 1 ) );
                A = + ( x_k - x_p ) * cos( phi ) + ( y_k - y_p ) * sin( phi );
                B = - ( x_k - x_p ) * sin( phi ) + ( y_k - y_p ) * cos( phi );
                C_x = y_k - y_p;
                C_y = x_k - x_p;
                D_x = sin( phi );
                D_y = cos( phi );
                L = norm( panel );
                world_influence_coeffs( 1 , i ) = - ( 1 / ( 2 * pi ) ) * ...
                    ( D_x / 2 * log( ( ( L + A ) ^ 2 + B ^ 2 ) / ( A ^ 2 + B ^ 2 ) ) + ...
                    ( C_x - A * D_x ) / B * ( atan( ( L + A ) / B ) - atan( A / B ) ) );
                world_influence_coeffs( 2 , i ) = + ( 1 / ( 2 * pi ) ) * ...
                    ( D_y / 2 * log( ( ( L + A ) ^ 2 + B ^ 2 ) / ( A ^ 2 + B ^ 2 ) ) + ...
                    ( C_y - A * D_y ) / B * ( atan( ( L + A ) / B ) - atan( A / B ) ) );
        end
    end
end

function [ norm_influence_coeffs , tan_influence_coeffs ] = gen_boundary_influence_coeffs( coordinates )

    n_panels = length( coordinates( : , 1 ) ) - 1;
    norm_influence_coeffs = zeros( n_panels );
    tan_influence_coeffs = zeros( n_panels );

    for i = 1 : n_panels

        x_p = ( coordinates( i , 1 ) + coordinates( i + 1 , 1 ) ) / 2;
        y_p = ( coordinates( i , 2 ) + coordinates( i + 1 , 2 ) ) / 2;
        panel = coordinates( i + 1 , : ) - coordinates( i , : );
        phi = atan2( panel( 2 ) , panel( 1 ) );
        rotation = [ cos( phi ) sin( phi ) ; -sin( phi ) cos( phi ) ];
        world_influence_coeffs = gen_world_influence_coeffs( coordinates , x_p , y_p );
        panel_influence_coeffs = rotation * world_influence_coeffs;
        tan_influence_coeffs( i , : ) = panel_influence_coeffs( 1 , : );
        norm_influence_coeffs( i , : ) = panel_influence_coeffs( 2 , : );

    end

    tan_influence_coeffs = tan_influence_coeffs + 0.5 * eye( n_panels );

end

function vortex_strengths = calculate_vortex_strengths( influence_coeffs , rel_velos )

    n_panels = length( rel_velos );
    influence_coeffs( end , : ) = zeros( [ 1 n_panels ] );
    influence_coeffs( end , 1 ) = 1;
    influence_coeffs( end , end ) = 1;
    rel_velos( end ) = 0;
    % influence_coeffs = influence_coeffs + 1e-5*eye(n_panels);
    vortex_strengths = influence_coeffs \ -rel_velos;

end

function kj_lift = kutta_joukowski_lift( coordinates , vortex_strengths )

    n_panels = length( coordinates( : , 1 ) ) - 1;
    kj_lift = 0;

    for i = 1 : n_panels
        panel_length = norm( coordinates( i + 1 , : ) - coordinates( i , : ) );
        kj_lift = kj_lift + 2 * panel_length * vortex_strengths( i );
    end

end

function pres_coeffs = pressure_coeffs( tan_influence_coeffs , vortex_strengths , tan_velos )

    total_tan_velos = tan_velos + ( tan_influence_coeffs - eye( length( tan_velos ) ) ) * vortex_strengths;
    pres_coeffs = 1 - total_tan_velos .^ 2;
    
end

function plot_pres_coeffs( coordinates , pres_coeffs )

   x_coords = ( coordinates( 1 : end - 1 , 1 ) + coordinates( 2 : end , 1 ) ) ./ 2;
   figure;
   hold on;
   grid on;
   plot( x_coords , pres_coeffs )
   set( gca , 'YDir' , 'reverse' )

end

function standard_lift = standard_lift( coordinates , pres_coeffs )

    [ ~ , min_x_index ] = min( coordinates( : , 1 ) );
    top_coords = coordinates( 1 : min_x_index , : );
    bot_coords = coordinates( min_x_index : end , : );
    top_dx_lengths = abs( top_coords( 2 : end , 1 ) - top_coords( 1 : end - 1 , 1 ) );
    bot_dx_lengths = abs( bot_coords( 2 : end , 1 ) - bot_coords( 1 : end - 1 , 1 ) );
    top_pres_coeffs = pres_coeffs( 1 : min_x_index - 1 );
    bot_pres_coeffs = pres_coeffs( min_x_index : end );
    standard_lift = sum( bot_dx_lengths .* bot_pres_coeffs ) - sum( top_dx_lengths .* top_pres_coeffs );
    
end


function streamline_plot( coordinates , vortex_strengths , attack_angle )

    figure;
    hold on;
    grid on;
    fill( coordinates( : , 1 ) , coordinates( : , 2 ) , 'k' );
    daspect( [ 1 1 1 ] );
    xlim( [ -0.75 , 1.75 ] );
    ylim( [ -0.5 , 0.5 ] );

    x_vals = -0.75 : 0.01 : 1.75;
    y_starts = -0.775 : 0.05 : 0.775;
    
    auxiliaries = [ [ vortex_strengths ; attack_angle ] coordinates ];

    for i = 1 : length( y_starts )
        [ x , y ] = runge_kutta_4( @streamline_derivative , x_vals , y_starts( i ) , auxiliaries );
        plot( x , y , 'b' );
    end

end

function streamline_slope = streamline_derivative( y , x , auxiliaries )
    coordinates = auxiliaries( : , 2 : 3 );
    vortex_strengths = auxiliaries( 1 : end - 1 , 1 );
    attack_angle = auxiliaries( end , 1 );
    influence_coeffs = gen_world_influence_coeffs( coordinates , x , y );
    velocity = influence_coeffs * vortex_strengths + [ cos( attack_angle * pi / 180 ) ; sin( attack_angle * pi / 180 ) ];
    streamline_slope = velocity( 2 ) / velocity( 1 );
end

function [ t , y ] = runge_kutta_4( system , times , initial_conditions , auxiliaries )
    t = times;
    y( 1 , : ) = initial_conditions;
    for i = 1 : ( length( t ) - 1 )
        h = t( i + 1 ) - t( i );
        K1 = h * system( y( i , : ) , t( i ) , auxiliaries );
        K2 = h * system( y( i , : ) + 0.5 * K1' , t( i ) + 0.5 * h , auxiliaries );
        K3 = h * system( y( i , : ) + 0.5 * K2' , t( i ) + 0.5 * h , auxiliaries );
        K4 = h * system( y( i , : ) + K3' , t( i ) + h , auxiliaries );
        changes = 1 / 6 * K1 + 1 / 3 * ( K2 + K3 ) + 1 / 6 * K4;
        y( i + 1 , : ) = y( i , : ) + changes';
    end
end

function pressure_contour_plot( coordinates , vortex_strengths , attack_angle )

    x_vals = -0.75: 0.01 : 1.75;
    y_vals = -0.75 : 0.01 : 0.75;
    [ X , Y ] = meshgrid( x_vals , y_vals );
    Z = zeros( [ length( y_vals ) length( x_vals ) ] );
    figure;
    hold on;
    xlim( [ -0.75 , 1.75 ] );
    ylim( [ -0.73 , 0.73 ] );
    daspect( [ 1 1 1 ] )
    freestream_velo = [ cos( attack_angle * pi / 180 ); sin( attack_angle * pi / 180 ) ];
    for i = 1 : length( y_vals )
        for j = 1 : length( x_vals )
            xy_influence_coeffs = gen_world_influence_coeffs( coordinates , X( i , j ) , Y( i , j ) );
            velocity = xy_influence_coeffs * vortex_strengths + freestream_velo;
            Z( i , j ) = 1 - norm( velocity ) ^ 2;
        end
    end

    contourf( X , Y , Z , 512 , 'EdgeColor' , 'none' );
    colormap( 'gray' );
    fill( coordinates( : , 1 ) , coordinates( : , 2 ) , 'k' );

end