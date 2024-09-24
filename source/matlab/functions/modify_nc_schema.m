function [m_ncs, vars_to_write] = modify_nc_schema(ncs, dat)
% A partially- or fully-populated MatLab file schema for a NetCDF-format file
%  is supplied as input, and a copy of it is then modified according to the
%  information in the supplied 'dat' structure.  Returns the modified file
%  schema and a cell array containing info on variables that may be written
%  later.

m_ncs = ncs;

% Set any attribute and dimension values, tabulate info on vars to write:
i_vw = 0;
vars_to_write = {};
f_fields = fieldnames(dat);
for i_f=1:length(f_fields)
   f_field = f_fields{i_f};
   if strcmp(f_field, 'global_atts')  % Global attribute(s)
      a_names = {m_ncs.Attributes.Name};
      a_fields = fieldnames(dat.(f_field));
      for i_a=1:length(a_fields)
         a_field = a_fields{i_a};
         a_idx = find(ismember(a_names, a_field));
         m_ncs.Attributes(a_idx).Value = dat.(f_field).(a_field);
      end
   else  % Group name
      g_names = {m_ncs.Groups.Name};
      g_idx = find(ismember(g_names, f_field));
      g_fields = fieldnames(dat.(f_field));
      for i_g=1:length(g_fields)
         g_field = g_fields{i_g};
         if strcmp(g_field, 'group_atts')  % Group attribute(s)
            a_names = {m_ncs.Groups(g_idx).Attributes.Name};
            a_fields = fieldnames(dat.(f_field).(g_field));
            for i_a=1:length(a_fields)
               a_field = a_fields{i_a};
               a_idx = find(ismember(a_names, a_field));
               m_ncs.Groups(g_idx).Attributes(a_idx).Value = ...
                                             dat.(f_field).(g_field).(a_field);
            end
         elseif strcmp(g_field, 'group_dims')  % Group dimension(s)
            d_names = {m_ncs.Groups(g_idx).Dimensions.Name};
            d_fields = fieldnames(dat.(f_field).(g_field));
            for i_d=1:length(d_fields)
               d_field = d_fields{i_d};
               d_idx = find(ismember(d_names, d_field));
               m_ncs.Groups(g_idx).Dimensions(d_idx).Length = ...
                                             dat.(f_field).(g_field).(d_field);
            end
            for i_v=1:length(m_ncs.Groups(g_idx).Variables)
               d_names = {m_ncs.Groups(g_idx).Variables(i_v).Dimensions.Name};
               for i_d=1:length(d_fields)
                  d_field = d_fields{i_d};
                  d_idx = find(ismember(d_names, d_field));
                  if ~isempty(d_idx)
                     m_ncs.Groups(g_idx).Variables(i_v).Dimensions(d_idx).Length = ...
                                             dat.(f_field).(g_field).(d_field);
                  end
               end
            end
         else  % Group variable(s)
            i_vw = i_vw+1;
            vars_to_write{i_vw} = {g_field, f_field, dat.(f_field).(g_field)};
         end
      end
   end
end

end
