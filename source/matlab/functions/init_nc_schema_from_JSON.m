function [ncs, jdat] = init_nc_schema_from_JSON(input_JSON_fpath, verbose)
% Returns a populated MatLab file schema for a NetCDF-format file, constructed
%  using information from a JSON-format data specification file.

% NumPy to MatLab datatype look-up table (structure):
NumPy_to_M_dtypes = struct;
NumPy_to_M_dtypes.float32 = 'single';
NumPy_to_M_dtypes.float64 = 'double';
NumPy_to_M_dtypes.int8 = 'int8';
NumPy_to_M_dtypes.int16 = 'int16';
NumPy_to_M_dtypes.int32 = 'int32';
NumPy_to_M_dtypes.int64 = 'int64';
NumPy_to_M_dtypes.uint8 = 'uint8';
NumPy_to_M_dtypes.uint16 = 'uint16';
NumPy_to_M_dtypes.uint32 = 'uint32';
NumPy_to_M_dtypes.uint64 = 'uint64';

%% Read in JSON-format filespec:
jdat = jsondecode(fileread(input_JSON_fpath));

%% Determine global attributes, and dimension names/types group attributes for
%%  every relevant group in the filespec:

global_Atts = {};
group_Atts = {};
i_global_Atts = 0;
i_group_Atts = 0;
chk_dims = struct;

group_names = fieldnames(jdat);
for i_g=1:length(group_names)
   g_name = group_names{i_g};

   if contains(g_name, 'JSON_COMMENTS')  % Contains only comments, so skip it
      if (verbose)
         sprintf('Irrelevant group %s -- skipping', g_name)
      end
      continue;

   elseif contains(g_name, '_Attributes')
      if strcmp(g_name, 'Global_Attributes')  % Contains global attribute(s)
         if length(jdat.(g_name)) > 0
            att_names = fieldnames(jdat.(g_name));
            for i_a=1:length(att_names)
               i_global_Atts = i_global_Atts+1;
               a_name = att_names{i_a};
               val = jdat.(g_name).(a_name).value;
               if isnumeric(val)
                  global_Atts{i_global_Atts} = {a_name, val};
               else
                  global_Atts{i_global_Atts} = {a_name, ...
                                              strrep(val, '!to_be_set!', ' ')};
               end
            end
         end
      else  % Contains group attribute(s)
         if length(jdat.(g_name)) > 0
            actual_g_name = strrep(g_name, '_Group_Attributes', '');
            att_names = fieldnames(jdat.(g_name));
            for i_a=1:length(att_names)
               i_group_Atts = i_group_Atts+1;
               a_name = att_names{i_a};
               val = jdat.(g_name).(a_name).value;
               if isnumeric(val)
                  group_Atts{i_group_Atts} = {actual_g_name, a_name, val};
               else
                  group_Atts{i_group_Atts} = {actual_g_name, a_name, ...
                                              strrep(val, '!to_be_set!', ' ')};
               end
            end
         end
      end

   else
      chk_dims.(g_name) = struct;
      var_names = fieldnames(jdat.(g_name));
      for i_v=1:length(var_names)
         v_name = var_names{i_v};
         dim_att = jdat.(g_name).(v_name).F_dimensions;
         for i_d=1:length(dim_att)
            tmp = dim_att{i_d};
            if contains(tmp, '!U!')
               d_name = strrep(tmp, '!U!', '');
               d_val = -1;
            else
               d_name = tmp;
               d_val = 0;
            end
            if ~isfield(chk_dims.(g_name), d_name)
               chk_dims.(g_name).(d_name) = {};
            end
            i = length(chk_dims.(g_name).(d_name))+1;
            chk_dims.(g_name).(d_name){i} = d_val;
         end
      end
   end
end

dims = struct;
group_names = fieldnames(chk_dims);
for i_g=1:length(group_names)
   g_name = group_names{i_g};
   dims.(g_name) = {};
   dim_names = fieldnames(chk_dims.(g_name));

   for i_d=1:length(dim_names)
      d_name = dim_names{i_d};

      if (length(chk_dims.(g_name).(d_name)) > 1)
         consistent = isequal(chk_dims.(g_name).(d_name){:});
      else
         consistent = true;
      end

      if (~consistent)
         fprintf(1, 'Type of dimension %s (group %s) was inconsistent.', ...
                                                               d_name, g_name);
      else
         i = length(dims.(g_name))+1;
         dims.(g_name){i} = {d_name, chk_dims.(g_name).(d_name){1}};
      end
   end
end

%% Build up the file schema structure:

% Begin with initializing the output file schema structure:
%not implemented here@ ncs.Filename = '';
ncs.Name = '/';
%not implemented here@ ncs.Dimensions = [struct];
%not implemented here@ ncs.Variables = [struct];
ncs.Groups = struct;
ncs.Format = 'netcdf4';

% Add the schema info for each global attribute:
if length(global_Atts) > 0
   ncs.Attributes = struct;
   for i_a=1:length(global_Atts)
      ncs.Attributes(i_a).Name = global_Atts{i_a}{1};
      ncs.Attributes(i_a).Value = global_Atts{i_a}{2};
   end
end

% Now add the schema info for each actual group:

i_rG = 0;
group_names = fieldnames(jdat);
for i_g=1:length(group_names)
   g_name = group_names{i_g};
   if contains(g_name, 'JSON_COMMENTS') || contains(g_name, '_Attributes')
      continue;  % Irrelevant group, skip it
   else
      i_rG = i_rG+1;
      ncs.Groups(i_rG).Name = g_name;
      ncs.Groups(i_rG).Dimensions = struct;
      ncs.Groups(i_rG).Variables = struct;

      % Group attributes:
      first_found = true;
      i_rA = 0;
      for i_a=1:length(group_Atts)
         a_ca = group_Atts{i_a};
         if strcmp(a_ca{1}, g_name)
            if first_found
               ncs.Groups(i_rG).Attributes = struct;
               first_found = false;
            end
            i_rA = i_rA+1;
            ncs.Groups(i_rG).Attributes(i_rA).Name = a_ca{2};
            ncs.Groups(i_rG).Attributes(i_rA).Value = a_ca{3};
         end
      end

      % Group dimensions:
      for i_d=1:length(dims.(g_name))
         d_ca = dims.(g_name){i_d};
         ncs.Groups(i_rG).Dimensions(i_d).Name = d_ca{1};
         val = d_ca{2};
         if val < 0
            ncs.Groups(i_rG).Dimensions(i_d).Length = 0;
            ncs.Groups(i_rG).Dimensions(i_d).Unlimited = true;
         else
            ncs.Groups(i_rG).Dimensions(i_d).Length = 1;
            ncs.Groups(i_rG).Dimensions(i_d).Unlimited = false;
         end
      end

      % Group variables:
      var_names = fieldnames(jdat.(g_name));
      for i_v=1:length(var_names)
         v_name = var_names{i_v};
         ncs.Groups(i_rG).Variables(i_v).Name = v_name;
         ncs.Groups(i_rG).Variables(i_v).Dimensions = struct;
         ncs.Groups(i_rG).Variables(i_v).Datatype = ...
                           NumPy_to_M_dtypes.(jdat.(g_name).(v_name).np_dtype);
         if isfield(jdat.(g_name).(v_name), 'fill_value')
            ncs.Groups(i_rG).Variables(i_v).FillValue = ...
                                             jdat.(g_name).(v_name).fill_value;
         end
         if isfield(jdat.(g_name).(v_name), 'compression')
            comp_a = jdat.(g_name).(v_name).compression;
            ncs.Groups(i_rG).Variables(i_v).DeflateLevel = comp_a{1};
            if comp_a{2} == true
               if verbose
                  fprintf(1, ['WARNING: NetCDF shuffle does not currently ' ...
                              'work properly with these MatLab tools']);
               end
%@               ncs.Groups(i_rG).Variables(i_v).Shuffle = logical(1);
%@            else
%@                ncs.Groups(i_rG).Variables(i_v).Shuffle = logical(0);
            end
         end

         % Variable dimensions:
         dim_att = jdat.(g_name).(v_name).F_dimensions;
         for i_d=1:length(dim_att)
            tmp = dim_att{i_d};
            if contains(tmp, '!U!')
               d_name = strrep(tmp, '!U!', '');
               ncs.Groups(i_rG).Variables(i_v).Dimensions(i_d).Length = 0;
               ncs.Groups(i_rG).Variables(i_v).Dimensions(i_d).Unlimited = true;
            else
               d_name = tmp;
               ncs.Groups(i_rG).Variables(i_v).Dimensions(i_d).Length = 1;
               ncs.Groups(i_rG).Variables(i_v).Dimensions(i_d).Unlimited = false;
            end
            ncs.Groups(i_rG).Variables(i_v).Dimensions(i_d).Name = d_name;
         end

         % Variable attributes:
         v_atts_to_ignore = {'np_dtype', 'C_dimensions', 'F_dimensions', ...
                             'fill_value', 'compression'};
         first_found = true;
         all_v_atts = fieldnames(jdat.(g_name).(v_name));
         i_u = 0; 
         for i_a=1:length(all_v_atts)
            a_name = all_v_atts{i_a};
            if ~any(strcmp(v_atts_to_ignore, a_name))
               if first_found
                  ncs.Groups(i_rG).Variables(i_v).Attributes = struct;
                  first_found = false;
               end
               i_u = i_u+1;
               ncs.Groups(i_rG).Variables(i_v).Attributes(i_u).Name = ...
                                                                        a_name;
               ncs.Groups(i_rG).Variables(i_v).Attributes(i_u).Value = ...
                                               jdat.(g_name).(v_name).(a_name); 
            end
         end
      end
   end
end

end
