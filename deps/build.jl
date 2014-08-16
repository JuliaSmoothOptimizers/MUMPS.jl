using BinDeps

@BinDeps.setup

@osx_only begin
  ChangeDirectory("src")
  run(`make`)
end
