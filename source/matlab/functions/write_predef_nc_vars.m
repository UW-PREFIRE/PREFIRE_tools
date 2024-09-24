function write_predef_nc_vars(nc_fpath, var_dat)
% Writes data to NetCDF-format file variables.  Assumes that the NetCDF-format
%  file has already been created and populated with metadata.  Variables' data
%  are passed in as a cell array -- each cell contains a cell array with info
%  about a variable to be written:
%    {1} = variable name (does not include any group name)
%    {2} = group name, if any
%    {3} = data to be written
% NOTE: Assumes start indices of 1 and a count vector equal to the
%        data array size().
 
ncid = netcdf.open(nc_fpath, 'WRITE');

for i_vw=1:length(var_dat)
   [v_name, g_name, v_data] = var_dat{i_vw}{:};

   if string(g_name) == ""
       gid = ncid;    
   else
       gid = netcdf.inqNcid(ncid, g_name);
   end
   varid = netcdf.inqVarID(gid, v_name);
   [o1, o2, dimids, o3] = netcdf.inqVar(gid, varid);
   n_dims = length(dimids);

   count = size(v_data);
   if (length(count) ~= n_dims)
      % Likely an "extraneous" dimension of length 1 in the 'v_data' array; use
      %  file spec dimension info instead:
      count = zeros(1, n_dims);
      for i_d=1:n_dims
         [o1, dimlen] = netcdf.inqDim(gid, dimids(i_d));
         count(i_d) = dimlen;
      end
   end
   start = zeros(1, n_dims);  % C-indexing for this API

   netcdf.putVar(gid, varid, start, count, v_data);
end

netcdf.close(ncid);

end
