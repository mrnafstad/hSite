<template>
  <v-app>
    <v-app-bar
      app
      color="dark"
      dark
    >
      Citybike-tracker
      <v-spacer />
      <v-btn @click="toggleAppInfo" >Info</v-btn>
    </v-app-bar>

    <v-main class="text-center">
      <GoogleMap />
    </v-main>
  </v-app>
</template>

<script>
import GoogleMap from './components/GoogleMap';
import { mapActions, mapGetters } from 'vuex'

export default {
  name: 'App',

  components: {
    GoogleMap,
  },
  computed: {
    ...mapGetters(['showAppInfo']),
  },
  methods: {
    ...mapActions(['toggleAppInfo'])
  },

  data: () => ({
    //
  }),
};
</script>
