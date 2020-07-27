%%    Minimal conductivity of graphene
%    Copyright (C) 2015 Peter Rakyta, Ph.D.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.
%
%> @addtogroup Examples
%> @{
%> @file MinimalConductivity.m
%> @brief Calculates the minimal conductivity of the graphene sheet as a function of the aspect ratio.
%> @param filenum The identification number of the filenema for the exported data (default is 1).
%> @Available
%> EQuUs v4.4 or later
%> @expectedresult
%> @image html MinimalConductivity.jpg
%> @image latex MinimalConductivity.jpg
%> The red line represents the theoretical \f$ 2/\pi\;\sigma_0 \f$ limit of the minimal conductivity with \f$ \sigma_0 = e^2/\hbar \f$.
%> @}
%
%> @brief Calculates the minimal conductivity of the graphene sheet as a function of the aspect ratio.
%> @param filenum The identification number of the filenema for the exported data (default is 1).
function ribbon_conductance( keep_sites )

if ~exist('filenum', 'var')
    filenum = 1;
end

filename = mfilename('fullpath');
[directory, fncname] = fileparts( filename );

    % filename containing the input XML    
    inputXML = 'Graphene_junction.xml';
    % Parsing the input file and creating data structures
    [Opt, param] = parseInput( fullfile( directory, inputXML ) );

    % filename containing the output XML
    outfilename = [fncname, '_',num2str( filenum )];
    % the output folder
    outputdir   = [];    

    % length of the junction (number of unit cells)
    %height = 10;
    % The Fermi energy
    EF = []; %In EV

    % class representing the two terminal ribbon structure
    cRibbon    = [];


    % creating 1D energy array
    % The minimal value of the widths
    width = 3;
    height=2;

    % creating 1D energy array (in units of eV)
    % minimum of the energy array
    Emin = 0.001;
    % maximum of the energy array
    Emax = 0.6;
    %> number of energy points
    Enum = 200;
    % the energy array
    Evec = Emin:(Emax-Emin)/(Enum+1):Emax;
    
    % creating output directories
    setOutputDir();

    % Calculate the transport through the magnetic barrier
    CalculateTransport();
    
%% CalculateTransport
%> @brief Calculates the conductance values
    function CalculateTransport()
        
        % determine the current energy value
        Energy = 1e-4; 

        %> creating the Ribbon class representing the twoterminal setup
        cRibbon = Ribbon('Opt', Opt, 'param', param, 'width', width, 'height', height );%,'leadmodel', @CustomLead);

        cRibbon.setEnergy( Energy )
        
        % Calculates the surface Green operator of the scattering region and perform a gauge transformation
        if keep_sites
            cRibbon.CalcFiniteGreensFunctionFromHamiltonian( 'gauge_trans', true ,'scatterPotential',@ScatterPot);      
        else
            cRibbon.CalcFiniteGreensFunctionFromHamiltonian( 'gauge_trans', true ,'scatterPotential',@ScatterPot);      
        end
        
        % Calculates the surface Green operator of the scattering region
        %cRibbon.CalcFiniteGreensFunction(); 

        % creating function handle for the Dyson Eq.       
        %Dysonfunc = @()cRibbon.CustomDysonFunc( 'constant_channels', false, 'SelfEnergy', selfEnergy , 'decimate', false);%, 'keep_sites','lead');

        [G,Ginv,junction_sites]= cRibbon.CustomDysonFunc( 'constant_channels', false, 'SelfEnergy', false , 'decimate', false,'UseHamiltonian',true);

        % Evaluate the Dyson Eq.
        %[G,junction_sites] = cRibbon.FL_handles.DysonEq( 'CustomDyson', Dysonfunc );
        %[a1,a2] = cRibbon.CalcSpectralFunction ( Energy , 'decimateDyson', false)
        
        cRibbon.CreateScatter();
        CreateH = cRibbon.CreateH;

        coordinates = CreateH.Read('coordinates');
        x = coordinates.x;
        y = coordinates.y;

        figure1 = figure('rend','painters','pos',[10 10 1200 1200]);       
        
        Ham = Ginv + eye(size(Ginv))*Energy;

        [V,D] = eig(Ham);
        
        Energies = diag(D);
        
        zero_energies_logical =  abs(Energies) < 2 ;
        
        zero_energies = Energies(zero_energies_logical);
        
        nr_zero_energies = sum(zero_energies_logical);
                
        zero_waveFsq = V(:,zero_energies_logical).*conj(V(:,zero_energies_logical));
        
        for j=1:nr_zero_energies               
            ax=subplot(ceil(sqrt(nr_zero_energies)),ceil(sqrt(nr_zero_energies)),j);
            plot(x,y,'x','color','black','MarkerSize',15)
            hold on
            scatter(x,y,zero_waveFsq(junction_sites.Scatter.site_indexes,j)*2000,'LineWidth',3,...
                'MarkerEdgeColor',[0 .5 .5],...
                'MarkerFaceColor',[0 .5 .5])%'.','color','red','MarkerSize',
            daspect([1 1 1]);
            xlim([min(x)-0.5,max(x)+0.5]);
            ylim([min(y)-0.5,max(y)+0.5]);
            title(['E = ',num2str(zero_energies(j))],'Interpreter','latex','FontSize',12);
        end    

        % exporting the calculated data
        print('-dpng', fullfile(outputdir,[outfilename,'W',num2str(width),'H',num2str(height),'keep_site',num2str(keep_sites),]) );
        disp(' Plotting. ')
        disp( ' ' )
        close(figure1);  
                 
    end
   function ret = ScatterPot( CreateH , Energy)
    
        coordinates = CreateH.Read('coordinates');
        x = coordinates.x;
        
        CreateH.Write('kulso_szabfokok',1:length(x));
        
        ret=zeros(1,length(x));
   end
%% setOutputDir
%> @brief Sets output directory.
    function setOutputDir()
        resultsdir = [pwd, filesep, 'results'];
        mkdir(resultsdir );
        outputdir = resultsdir;        
    end

    function lead_tmp=CustomLead(idx, E, varargin)

        p_inner = inputParser;
        
        p_inner.addParameter('createCore', 0);
        p_inner.addParameter('Just_Create_Hamiltonians', 0);
        p_inner.addParameter('shiftLead', 0);
        p_inner.addParameter('coordinates_shift', 0);
        p_inner.addParameter('transversepotential', []);
        p_inner.addParameter('Lead',[]);
        p_inner.addParameter('gauge_field', [] );% gauge field for performing gauge transformation
        p_inner.addParameter('SelfEnergy',false); % set true to calculate the self energy of the semi-infinite lead
        p_inner.addParameter('SurfaceGreensFunction', true );% set true to calculate the surface Greens function of the semi-infinite lead
        p_inner.addParameter('q', [] )
        
        p_inner.parse(varargin{:});
        
        createCore            = p_inner.Results.createCore;
        Just_Create_Hamiltonians = p_inner.Results.Just_Create_Hamiltonians;
        shiftLead             = p_inner.Results.shiftLead;
        coordinates_shift     = p_inner.Results.coordinates_shift;
        transversepotential   = p_inner.Results.transversepotential;
        Lead                  = [];%p_inner.Results.Surface_tmp;
        gauge_field           = p_inner.Results.gauge_field; % gauge field for performing gauge transformation
        SelfEnergy            = p_inner.Results.SelfEnergy;
        SurfaceGreensFunction = p_inner.Results.SurfaceGreensFunction;
        q                     = p_inner.Results.q;

        % setting th ephysical values of the square lattice
        
        Opt_lead = Opt;
        param_lead = param;
                
        Opt_lead.Lattice_Type = 'Graphene2';
        param_lead.Leads{idx} = param_Graphene_Lead();

        [Opt_lead, param_lead] = ValidateStructures( Opt_lead, param_lead );    
        
        param_lead.Leads{idx}.M=width;
        param_lead.Leads{idx}.End_Type='Z';        
        
        param_lead.Leads{idx}.Lead_Orientation= -(-1)^idx;
        
        FL_handles = Transport_Interface(E, Opt_lead, param_lead);
        
        lead_tmp = FL_handles.SurfaceGreenFunctionCalculator(idx, 'createCore', createCore, ...
                            'Just_Create_Hamiltonians', Just_Create_Hamiltonians, ...
                            'shiftLead', shiftLead, ...
                            'coordinates_shift', coordinates_shift, ...
                            'transversepotential', transversepotential, ...
                            'gauge_field', gauge_field, ...
                            'SelfEnergy', SelfEnergy, ...
                            'SurfaceGreensFunction', SurfaceGreensFunction, ...
                            'q', q);   
    end
   
