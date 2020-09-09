<template>
  <v-container>
     <v-overlay :value="showAppInfo">
        <v-card raised light>
            This is a tracking app for Oslo city bikes
            <br />
        <v-btn icon @click="toggleAppInfo">
          <v-icon>mdi-check</v-icon>
        </v-btn>
      </v-card>
    </v-overlay>

    <v-overlay :value="showStationInfo">
      <v-card raised light>
        <v-list-item two-line>
          <v-list-item-content>
        <v-list-item-title>{{ currentInfo.name }}</v-list-item-title>
        <v-list-item-subtitle>
        Ledig: {{ currentInfo.available }}/{{ currentInfo.capacity }}
        </v-list-item-subtitle>
        <v-spacer></v-spacer>
        <v-btn icon @click="toggleStationInfo">
          <v-icon>mdi-check</v-icon>
        </v-btn>
        </v-list-item-content>
        </v-list-item>
      </v-card>
    </v-overlay>
    <v-overlay :value="user.show">
      <v-card raised light>
        <v-list-item>
          
        <v-list-item-title>{{ currentInfo.name }}</v-list-item-title>
        <v-spacer></v-spacer>
        <v-btn icon @click="user.show = !user.show">
          <v-icon>mdi-check</v-icon>
        </v-btn>
        </v-list-item>
      </v-card>
    </v-overlay>
    <gmap-map
      :center="center"
      :zoom="12"
      style="width:100%;  height: 90vh;"
    >
      <gmap-marker
        :key="index"
        v-for="(m, index) in stations"
        :position="m.pos"
        @click="toggleStation(m)"
        :clickable="true"
      >
      </gmap-marker>
    <gmap-marker 
      :position="center" 
      :clickable="true"
      @click="user.show = !user.show"
      >
    </gmap-marker>
    </gmap-map>
  </v-container>
</template>

<script>
import { mapActions, mapGetters } from 'vuex'
export default {
  name: "GoogleMap",
  data() {
    return {
        // default to Montreal to keep it simple
        // change this to whatever makes sense
        center: { lat: 45.508, lng: -73.587 },
        currentPlace: null,
        showInfo: true,
        currentInfo: {
          name: 'no name',
          capacity: 'no data'
        },
        user: {
          show: false,
          info: 'Users position'
        }
    };
  },

  mounted() {
    this.geolocate();
  },

  methods: {
      ...mapActions(['getStations', 'getStationData', 'toggleAppInfo', 'toggleStationInfo']),
    // receives a place object via the autocomplete component
    setPlace(place) {
      this.currentPlace = place;
    },
    geolocate: function() {
      navigator.geolocation.getCurrentPosition(position => {
        this.center = {
          lat: position.coords.latitude,
          lng: position.coords.longitude
        };
      });
      //console.log(this.center)
    },
    toggleStation(station) {
      this.toggleStationInfo()
      this.currentInfo = {
        name: station.name,
        capacity: station.capacity,
        available: station.available,
        docksFree: station.docksAvailable
      }
    }
  },
  computed: {
      ...mapGetters(['stations', 'showAppInfo', 'showStationInfo'])
  },
    async created() {
      await this.getStations()
      //console.log(this.stations)
      await this.getStationData()
      this.timer = setInterval(this.getStationData, 60)
  }
};
</script>