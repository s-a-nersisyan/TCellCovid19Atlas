$(window).on("load", function() {
  $('body').scrollspy({ target: '#navbar-protein', offset: 60 });

  $(window).on('activate.bs.scrollspy', function () {
    history.replaceState({}, "", $('.nav-item .active').attr("href"));
  });
});
